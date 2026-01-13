/********************************************************************************
 * SPDX-License-Identifier: EPL-2.0
 ********************************************************************************/

#include "planning/speed_profile.hpp"

#include <cmath>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <OsqpEigen/OsqpEigen.h>

#include "adore_map/map_point.hpp"
#include "adore_math/curvature.hpp"
#include "adore_math/distance.h"
#include "adore_math/point.h"
#include "adore_math/geometry/projection.hpp"
#include "adore_math/geometry/interval.hpp"

namespace adore::planner
{
namespace
{

constexpr double k_inf = 1e30;

// =============================================================================
// Internal Types (formerly st_types.hpp)
// =============================================================================

enum class StBoundaryType : uint8_t
{
  Unknown = 0,
  DynamicObstacle,
  StaticObstacle,
};

struct StBoundarySample
{
  double t       = 0.0;
  double s_lower = 0.0;
  double s_upper = 0.0;
};

struct StBoundary
{
  int64_t                       obstacle_id = -1;
  StBoundaryType                type        = StBoundaryType::Unknown;
  std::vector<StBoundarySample> samples;
};

// =============================================================================
// ST Graph Builder Logic
// =============================================================================

StBoundaryType classify_boundary_type( const dynamics::TrafficParticipant& p )
{
  const double v = std::hypot( p.state.vx, p.state.vy );
  // Simple heuristic: if slow, static.
  if( v < 0.1 )
  {
    return StBoundaryType::StaticObstacle;
  }
  return StBoundaryType::DynamicObstacle;
}

std::optional<StBoundary>
build_boundary_for_participant( const DrivableArea& area, const dynamics::TrafficParticipant& p,
                                 const dynamics::PhysicalVehicleParameters& vehicle_params,
                                 const SpeedProfileConfig& cfg,
                                 double ego_s_seed, double t0_abs )
{
  StBoundary boundary;
  boundary.obstacle_id = p.id;
  boundary.type        = classify_boundary_type( p );

  // Sampling
  const double t_horizon = cfg.total_time;
  const double dt        = cfg.obstacle_dt; 
  const int    n_samples = static_cast<int>( std::ceil( t_horizon / dt ) );

  bool any_overlap = false;

  for( int i = 0; i <= n_samples; ++i )
  {
    const double t_rel = i * dt;

    const double v = std::hypot( p.state.vx, p.state.vy );
    const double yaw = p.state.yaw_angle;

    double px = p.state.x + v * std::cos( yaw ) * t_rel;
    double py = p.state.y + v * std::sin( yaw ) * t_rel;
    
    const auto& corners0 = p.get_corners().points;
    if (corners0.empty()) continue;
    
    double s_min = std::numeric_limits<double>::infinity();
    double s_max = -std::numeric_limits<double>::infinity();
    double l_min = std::numeric_limits<double>::infinity();
    double l_max = -std::numeric_limits<double>::infinity();

    for (const auto& c0 : corners0) {
        double dx = c0.x - p.state.x;
        double dy = c0.y - p.state.y;
        double cx = px + dx;
        double cy = py + dy;

        auto sl = adore::math::project_to_path( area.reference_line, adore::math::Point2d{cx, cy}, std::optional<double>(ego_s_seed), cfg.projection_window );
        if (sl) {
            s_min = std::min(s_min, sl->s);
            s_max = std::max(s_max, sl->s);
            l_min = std::min(l_min, sl->l);
            l_max = std::max(l_max, sl->l);
        }
    }

    if (s_min > s_max) continue; 

    // Inflate
    s_min -= cfg.obstacle_longitudinal_buffer;
    s_max += cfg.obstacle_longitudinal_buffer;
    l_min -= cfg.obstacle_lateral_buffer;
    l_max += cfg.obstacle_lateral_buffer;

    double s_avg = 0.5 * (s_min + s_max);
    auto it_left = area.left_boundary.lower_bound(s_avg);
    auto it_right = area.right_boundary.lower_bound(s_avg);
    
    if (it_left == area.left_boundary.end() || it_right == area.right_boundary.end()) continue;

    auto sl_left = adore::math::project_to_path(area.reference_line, adore::math::Point2d{it_left->second.x, it_left->second.y}, std::optional<double>(s_avg), 10.0);
    auto sl_right = adore::math::project_to_path(area.reference_line, adore::math::Point2d{it_right->second.x, it_right->second.y}, std::optional<double>(s_avg), 10.0);
    
    if (!sl_left || !sl_right) continue;

    // Check if obstacle overlaps with the vehicle footprint
    const double safety_margin = ( vehicle_params.body_width * 0.5 ) + cfg.obstacle_lateral_buffer;
    
    const double overlap_min = std::max( l_min, -safety_margin );
    const double overlap_max = std::min( l_max, safety_margin );

    if (overlap_min < overlap_max) {
        StBoundarySample sample;
        sample.t = t_rel;
        sample.s_lower = s_min;
        sample.s_upper = s_max;
        boundary.samples.push_back(sample);
        any_overlap = true;
    }
  }

  if (any_overlap) {
      return boundary;
  }
  return std::nullopt;
}

std::vector<std::vector<adore::math::Interval<double>>>
build_forbidden_intervals_per_t( const std::vector<double>& t_samples, const std::vector<StBoundary>& boundaries )
{
  std::vector<std::vector<adore::math::Interval<double>>> result( t_samples.size() );

  for( size_t k = 0; k < t_samples.size(); ++k )
  {
    double t = t_samples[k];
    for( const auto& b : boundaries )
    {
      if (b.samples.empty()) continue;
      
      auto it = std::lower_bound(b.samples.begin(), b.samples.end(), t, [](const StBoundarySample& s, double val){
          return s.t < val;
      });
      
      double s_lo = 0, s_hi = 0;
      bool valid = false;

      if (it == b.samples.begin()) {
          if (std::abs(it->t - t) < 0.1) {
             s_lo = it->s_lower; s_hi = it->s_upper; valid = true;
          }
      } else if (it == b.samples.end()) {
          auto prev = std::prev(it);
          if (std::abs(prev->t - t) < 0.1) {
             s_lo = prev->s_lower; s_hi = prev->s_upper; valid = true;
          }
      } else {
          auto prev = std::prev(it);
          double r = (t - prev->t) / (it->t - prev->t);
          s_lo = prev->s_lower + r * (it->s_lower - prev->s_lower);
          s_hi = prev->s_upper + r * (it->s_upper - prev->s_upper);
          valid = true;
      }

      if (valid) {
          result[k].push_back( adore::math::Interval<double>( s_lo, s_hi ) );
      }
    }
    adore::math::merge_intervals_in_place( result[k] );
  }
  return result;
}

// =============================================================================
// Curvature + legal speed -> reference SpeedProfile along s
// =============================================================================

// Conservative vmax lookup from speed_limit_profile_s (min of bracketing v_limit samples)
double
v_limit_at_s_min_bracket( const SpeedProfile& speed_limit_profile_s, double s_query, double default_v )
{
  if( speed_limit_profile_s.empty() )
  {
    return default_v;
  }

  if( s_query <= speed_limit_profile_s.front().s )
  {
    return speed_limit_profile_s.front().v_limit;
  }
  if( s_query >= speed_limit_profile_s.back().s )
  {
    return speed_limit_profile_s.back().v_limit;
  }

  const auto it = std::lower_bound( speed_limit_profile_s.begin(), speed_limit_profile_s.end(), s_query,
                                    []( const SpeedProfilePoint& p, double s_val ) { return p.s < s_val; } );

  const size_t i1 = static_cast<size_t>( std::distance( speed_limit_profile_s.begin(), it ) );
  const size_t i0 = i1 > 0 ? i1 - 1 : 0;

  return std::min( speed_limit_profile_s[i0].v_limit, speed_limit_profile_s[i1].v_limit );
}

double
kappa_at_s_nearest( const SpeedProfile& speed_limit_profile_s, double s_query )
{
  if( speed_limit_profile_s.empty() )
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const auto it = std::lower_bound( speed_limit_profile_s.begin(), speed_limit_profile_s.end(), s_query,
                                    []( const SpeedProfilePoint& p, double s_val ) { return p.s < s_val; } );

  if( it == speed_limit_profile_s.begin() )
  {
    return it->kappa;
  }
  if( it == speed_limit_profile_s.end() )
  {
    return speed_limit_profile_s.back().kappa;
  }

  const auto& p1 = *it;
  const auto& p0 = *( it - 1 );
  return ( std::abs( p1.s - s_query ) < std::abs( s_query - p0.s ) ) ? p1.kappa : p0.kappa;
}

// =============================================================================
// DP Optimizer Logic
// =============================================================================

struct DpCell
{
  double cost = k_inf;
  int    prev = -1;
  double v    = 0.0;
  double a    = 0.0;
  bool   ok   = false;
};

double obstacle_proximity_cost( const SpeedProfileConfig& cfg, double s, double t, const std::vector<adore::math::Interval<double>>& forbidden )
{
    double min_dist = k_inf;
    for (const auto& inv : forbidden) {
        double dist = 0.0;
        if (s < inv.lo) dist = inv.lo - s;
        else if (s > inv.hi) dist = s - inv.hi;
        else dist = 0.0;
        min_dist = std::min(min_dist, dist);
    }
    
    if (min_dist < 1e-3) return cfg.w_dp_obstacle_proximity * 100.0;
    
    return cfg.w_dp_obstacle_proximity * std::exp( -min_dist );
}

SpeedProfile
solve_dp( const SpeedProfileConfig& cfg,
          const dynamics::ComfortSettings& comfort,
          double ego_s, double ego_v, double ego_a,
          const std::vector<double>& t_samples,
          const std::vector<std::vector<adore::math::Interval<double>>>& forbidden_intervals,
          const SpeedProfile& v_limit_profile )
{
  const size_t N_t = t_samples.size();
  std::cerr << "solve_dp: ego_s=" << ego_s << ", ego_v=" << ego_v << ", a_max=" << comfort.max_acceleration << ", a_min=" << comfort.min_acceleration << ", v_max=" << comfort.max_speed << std::endl;
  if( N_t < 2 )
    return {};

  const double s_start = ego_s;
  const double s_end = ego_s + cfg.s_horizon;
  const double ds = cfg.ds_dp;
  const int N_s = static_cast<int>(std::ceil((s_end - s_start) / ds));
  
  std::vector<std::vector<DpCell>> cost_table( N_t, std::vector<DpCell>( N_s + 1 ) );

  cost_table[0][0].cost = 0.0;
  cost_table[0][0].v = ego_v;
  cost_table[0][0].a = ego_a;
  cost_table[0][0].ok = true;

  for( size_t k = 0; k < N_t - 1; ++k )
  {
    double dt = t_samples[k+1] - t_samples[k];
    
    for( int i = 0; i <= N_s; ++i )
    {
      if( !cost_table[k][i].ok ) continue;
      
      const double s_curr = s_start + i * ds;
      const double v_curr = cost_table[k][i].v;
      
      double s_next_min = s_curr + std::max(0.0, v_curr * dt + 0.5 * comfort.min_acceleration * dt * dt);
      double s_next_max = s_curr + std::max(0.0, v_curr * dt + 0.5 * comfort.max_acceleration * dt * dt);
      
      int j_min = static_cast<int>((s_next_min - s_start) / ds);
      int j_max = static_cast<int>((s_next_max - s_start) / ds) + 1;
      
      j_min = std::max(i, j_min);
      j_max = std::min(N_s, j_max);
      
      for( int j = j_min; j <= j_max; ++j )
      {
          double s_next = s_start + j * ds;
          double dist = s_next - s_curr;
          double v_next = dist / dt * 2 - v_curr;
          double a = (v_next - v_curr) / dt;
          
          const double local_v_limit = v_limit_at_s_min_bracket( v_limit_profile, s_next, comfort.max_speed );

          if (v_next < cfg.v_min - 1e-3 || v_next > comfort.max_speed + 1e-3) continue;
          if (v_next > local_v_limit + 1e-3) continue;
          if (a < comfort.min_acceleration - 1e-3 || a > comfort.max_acceleration + 1e-3) continue;
          
          bool collision = false;
          if (k + 1 < forbidden_intervals.size()) {
              for( const auto& inv : forbidden_intervals[k+1] ) {
                  if (inv.contains(s_next)) { collision = true; break; }
              }
          }
          if (collision) continue;

          double edge_cost = 0.0;
          
          // Target speed: if a specific cruise speed is set, use it (clipped by local limit); 
          // otherwise target the local limit (which includes map limits and curvature caps).
          double target_v = (cfg.v_cruise > 0.0) ? std::min(cfg.v_cruise, local_v_limit) : local_v_limit;

          edge_cost += cfg.w_dp_speed * std::abs(v_next - target_v);
          edge_cost += cfg.w_dp_accel * a * a;
          edge_cost += cfg.w_dp_jerk * std::abs(a - cost_table[k][i].a) / dt; 
          edge_cost += -cfg.w_dp_progress * dist;
          
          if (k+1 < forbidden_intervals.size())
             edge_cost += obstacle_proximity_cost(cfg, s_next, t_samples[k+1], forbidden_intervals[k+1]);

          double total = cost_table[k][i].cost + edge_cost;
          
          if (total < cost_table[k+1][j].cost) {
              cost_table[k+1][j].cost = total;
              cost_table[k+1][j].prev = i;
              cost_table[k+1][j].v = v_next;
              cost_table[k+1][j].a = a;
              cost_table[k+1][j].ok = true;
          }
      }
    }
  }

  int best_end_idx = -1;
  double min_c = k_inf;
  for( int i = 0; i <= N_s; ++i ) {
      if (cost_table[N_t-1][i].ok && cost_table[N_t-1][i].cost < min_c) {
          min_c = cost_table[N_t-1][i].cost;
          best_end_idx = i;
      }
  }

  if (best_end_idx == -1) {
      std::cerr << "solve_dp FAILED: no valid end state found" << std::endl;
      return {};
  }

  SpeedProfile result;
  result.resize(N_t);
  
  int curr = best_end_idx;
  for (int k = N_t - 1; k >= 0; --k) {
      result[k].t = t_samples[k];
      result[k].s = s_start + curr * ds;
      result[k].v = cost_table[k][curr].v;
      result[k].a = cost_table[k][curr].a;
      
      curr = cost_table[k][curr].prev;
  }
  
  return result;
}

// =============================================================================
// QP Smoother Logic
// =============================================================================

SpeedProfile
smooth_qp( const SpeedProfileConfig& cfg, double ego_s, double ego_v, double ego_a,
         const SpeedProfile& dp_profile,
         const std::vector<std::vector<adore::math::Interval<double>>>& forbidden_intervals )
{
  if (dp_profile.size() < 2) return dp_profile;

  // Placeholder for QP smoothing implementation
  // Currently returning DP profile directly as it provides a feasible solution.
  // Full QP implementation would refine smoothness and strictly enforce constraints.
  
  return dp_profile; 
}

// =============================================================================
// Small helpers for main pipeline
// =============================================================================

std::vector<double>
build_time_samples( double total_time, double dt )
{
  const int           n = std::max( 1, static_cast<int>( std::ceil( total_time / dt ) ) );
  std::vector<double> t( static_cast<size_t>( n ) + 1u, 0.0 );
  for( int i = 0; i <= n; ++i )
  {
    t[static_cast<size_t>( i )] = static_cast<double>( i ) * dt;
  }
  return t;
}

std::optional<double>
ego_s_from_area( const DrivableArea& area, const dynamics::VehicleStateDynamic& ego, const SpeedProfileConfig& cfg )
{
  adore::math::Point2d p;
  p.x = ego.x;
  p.y = ego.y;

  const auto sl = adore::math::project_to_path( area.reference_line, p, std::optional<double>(std::nullopt), cfg.projection_window );
  if( !sl )
  {
    return std::nullopt;
  }
  return sl->s;
}

// =============================================================================
// Small helpers for main pipeline
// =============================================================================

void
speed_forward_backward_pass( std::vector<SpeedProfilePoint>& profile, double ds, double a_max, double a_min, double v_min, double v_max,
                             double v0_hard )
{
  if( profile.size() < 2 || ds <= 1e-6 )
    return;

  const double a_up   = std::max( 0.0, a_max );
  const double a_down = std::max( 0.0, -a_min );

  auto clamp_v = [&]( double v ) { return std::clamp( v, std::max( 0.0, v_min ), v_max ); };

  // --- 1) Start from caps = profile[i].v_limit and make them brake-feasible (backward tighten)
  std::vector<double> v_cap( profile.size(), 0.0 );
  for( size_t i = 0; i < profile.size(); ++i )
  {
    v_cap[i] = clamp_v( profile[i].v_limit );
  }

  if( a_down > 0.0 )
  {
    for( size_t i = profile.size() - 1; i > 0; --i )
    {
      const double v_i       = std::max( 0.0, v_cap[i] );
      const double vmax_prev = std::sqrt( std::max( 0.0, v_i * v_i + 2.0 * a_down * ds ) );
      v_cap[i - 1]           = std::min( v_cap[i - 1], clamp_v( vmax_prev ) );
    }
  }

  // --- 2) Forward propagate actual profile with hard v0 and "max-brake fallback" when cap is unreachable
  std::vector<double> v( profile.size(), 0.0 );
  v[0] = clamp_v( v0_hard ); // HARD constraint

  for( size_t i = 0; i + 1 < v.size(); ++i )
  {
    const double v_i = std::max( 0.0, v[i] );

    const double v_max_acc   = ( a_up > 0.0 ) ? std::sqrt( std::max( 0.0, v_i * v_i + 2.0 * a_up * ds ) ) : v_i;
    const double v_min_brake = ( a_down > 0.0 ) ? std::sqrt( std::max( 0.0, v_i * v_i - 2.0 * a_down * ds ) ) : v_i;

    const double v_hi = std::min( v_cap[i + 1], clamp_v( v_max_acc ) );

    // If reachable under decel limit, respect cap; otherwise brake as hard as possible.
    v[i + 1] = ( v_hi >= v_min_brake ) ? v_hi : clamp_v( v_min_brake );
  }

  // write back
  for( size_t i = 0; i < profile.size(); ++i )
  {
    profile[i].v_limit = v_cap[i]; // keep the tightened cap if you want for debugging
    profile[i].v       = v[i];
  }

  profile.front().t = 0.0;
  for( size_t i = 0; i + 1 < profile.size(); ++i )
  {
    const double v_i   = std::max( 0.0, clamp_v( profile[i].v ) );
    const double v_ip1 = std::max( 0.0, clamp_v( profile[i + 1].v ) );
    const double dt    = ( std::abs( v_i + v_ip1 ) < 1e-6 ) ? 0.0 : ( 2.0 * ds ) / ( v_i + v_ip1 );
    profile[i + 1].t   = profile[i].t + dt;
  }
}

SpeedProfile
build_speed_limit_profile_s( const DrivableArea& area,
                             const dynamics::PhysicalVehicleParameters& vehicle_params,
                             const dynamics::ComfortSettings& comfort_settings,
                             const SpeedProfileConfig& cfg, double v0 )
{
  std::cerr << "build_speed_limit_profile_s: v0=" << v0 << ", comfort_v_max=" << comfort_settings.max_speed << std::endl;
  SpeedProfile ref_profile;

  const auto& ref = area.reference_line;
  if( ref.size() < 2 )
  {
    return ref_profile;
  }

  const double s_min = ref.begin()->first;
  const double s_max = ref.rbegin()->first;

  const double ds     = std::max( 0.2, cfg.speed_limit_ds );
  const double span_s = std::max( 0.0, s_max - s_min );
  const size_t n      = std::max<size_t>( 1u, static_cast<size_t>( std::ceil( span_s / ds ) ) );

  ref_profile.reserve( n + 1u );

  auto has_speed = []( const adore::map::MapPoint& mp ) {
    return mp.max_speed.has_value() && std::isfinite( *mp.max_speed ) && *mp.max_speed >= 0.0;
  };

  // Curvature: maintain "lower_bound(s)" iterator by only advancing (s is monotone).
  auto it_center = ref.begin();

  // Legal speed: nearest previous/next speed-marked points, advanced linearly.
  auto it_prev_speed = ref.end();
  auto it_next_speed = ref.begin();
  while( it_next_speed != ref.end() && !has_speed( it_next_speed->second ) )
  {
    ++it_next_speed;
  }

  auto legal_speed_from_bracketing = [&]( double default_v ) -> double {
    if( it_prev_speed == ref.end() && it_next_speed == ref.end() )
    {
      return default_v;
    }
    if( it_prev_speed == ref.end() )
    {
      return *it_next_speed->second.max_speed;
    }
    if( it_next_speed == ref.end() )
    {
      return *it_prev_speed->second.max_speed;
    }
    return std::min( *it_prev_speed->second.max_speed, *it_next_speed->second.max_speed );
  };

  auto curvature_at_s_monotonic = [&]( double s ) -> double {
    if( ref.size() < 3 )
      return 0.0;

    // advance center monotonically
    while( it_center != ref.end() && it_center->first < s )
      ++it_center;

    if( it_center == ref.end() )
      it_center = std::prev( ref.end() );

    // choose it1 as a valid interior point
    auto it1 = it_center;
    if( it1 == ref.begin() )
      it1 = std::next( it1 );
    if( std::next( it1 ) == ref.end() )
      it1 = std::prev( it1 );

    const double s1 = it1->first;

    // configurable stencil half-width in s (meters)
    constexpr double ds_kappa = 2.0;

    // find it0 at least ds_kappa behind
    auto it0 = ref.lower_bound( s1 - ds_kappa );
    // lower_bound gives first >= (s1-ds); we want <=, so step back if possible
    if( it0 == ref.end() )
      it0 = std::prev( ref.end() );
    if( it0->first > s1 - ds_kappa && it0 != ref.begin() )
      it0 = std::prev( it0 );

    // find it2 at least ds_kappa ahead
    auto it2 = ref.lower_bound( s1 + ds_kappa );
    if( it2 == ref.end() )
      it2 = std::prev( ref.end() );

    // ensure ordering and distinct iterators
    if( it0 == it1 && it0 != ref.begin() )
      it0 = std::prev( it0 );
    if( it2 == it1 && std::next( it2 ) != ref.end() )
      it2 = std::next( it2 );

    // final safety: need three distinct points in increasing order
    if( it0 == it1 || it1 == it2 || it0 == it2 )
      return 0.0;
    if( it0->first >= it1->first || it1->first >= it2->first )
      return 0.0;

    // extra guard against too-tiny spacing (even after ds_kappa, duplicates may exist) using math::distance_2d
    const double d01 = adore::math::distance_2d( it0->second, it1->second );
    const double d12 = adore::math::distance_2d( it1->second, it2->second );
    if( d01 < 1.0 || d12 < 1.0 )
      return 0.0;

    const double kappa = adore::math::compute_curvate( it0->second, it1->second, it2->second );
    
    // Derive max_curvature from vehicle parameters if it exceeds physical limits
    // kappa = 1/R, R = wheelbase / tan(delta_max)
    const double max_steer = std::max(std::abs(vehicle_params.steering_angle_max), std::abs(vehicle_params.steering_angle_min));
    const double R_min = vehicle_params.wheelbase / std::max(0.01, std::tan(max_steer));
    const double physical_max_curvature = 1.0 / std::max(0.1, R_min);
    const double effective_max_curvature = std::min(cfg.max_curvature, physical_max_curvature);

    if( kappa > effective_max_curvature * 3.0 )
    {
      return 0.0;
    }
    return std::isfinite( kappa ) ? std::abs( kappa ) : 0.0;
  };


  for( int i = 0; i <= n; ++i )
  {
    const double s_raw = s_min + static_cast<double>( i ) * ds;
    const double s     = std::clamp( s_raw, s_min, s_max );

    while( it_next_speed != ref.end() && it_next_speed->first <= s )
    {
      it_prev_speed = it_next_speed;
      ++it_next_speed;
      while( it_next_speed != ref.end() && !has_speed( it_next_speed->second ) )
      {
        ++it_next_speed;
      }
    }

    double kappa = curvature_at_s_monotonic( s );
    kappa        = std::min( kappa, cfg.max_curvature );

    double v_legal = legal_speed_from_bracketing( cfg.default_legal_speed );
    v_legal        = std::max( 0.0, v_legal ) * comfort_settings.speed_fraction_of_limit;
    v_legal        = std::clamp( v_legal, std::max( 0.0, cfg.v_min ), comfort_settings.max_speed );

    SpeedProfilePoint p;
    p.s       = s;
    p.t       = 0.0;
    p.kappa   = kappa;   // raw
    p.v_limit = v_legal; // legal for now
    if (v_legal < 0.1) std::cerr << " WARNING: v_legal=" << v_legal << " at s=" << s << std::endl;
    ref_profile.push_back( p );
  }

  int i = 0;
  for( auto& p : ref_profile )
  {
    if( p.kappa < cfg.curvature_eps )
    {
      p.kappa = 0.0;
      continue;
    }
    const double v_curv = std::sqrt( std::max( 0.0, comfort_settings.max_lateral_acceleration ) / p.kappa );
    const double v      = std::min( p.v_limit, v_curv );
    if (i % 10 == 0) std::cerr << " i=" << i << " s=" << p.s << " v_legal=" << p.v_limit << " v_curv=" << v_curv << " final_v_limit=" << v << std::endl;
    p.v_limit           = v;
    ++i;
  }
  const double s_end   = ref_profile.back().s;
  const double s_start = ref_profile.front().s;
  const double span    = s_end - s_start;
  if( span < cfg.s_horizon )
  {
    const size_t n_tail = std::min<size_t>( 5u, ref_profile.size() );
    for( size_t k = 0; k < n_tail; ++k )
    {
      const size_t i         = ref_profile.size() - 1u - k;
      ref_profile[i].v_limit = 0.0;
    }
  }


  return ref_profile;
}

} // namespace


// =============================================================================
// Public API Implementation
// =============================================================================

SpeedProfile plan_speed_profile( const DrivableArea& area, const dynamics::TrafficParticipantSet& participants,
                                 const dynamics::VehicleStateDynamic& ego_state,
                                 const dynamics::PhysicalVehicleParameters& vehicle_params,
                                 const dynamics::ComfortSettings& comfort_settings,
                                 const SpeedProfileConfig& config )
{
  if( area.empty() )
  {
    return {};
  }

  // 1. Ego S
  auto s_ego_opt = ego_s_from_area( area, ego_state, config );
  if( !s_ego_opt )
  {
    return {};
  }
  double ego_s = *s_ego_opt;
  double ego_v = std::hypot( ego_state.vx, ego_state.vy );

  // 2. Speed Limits
  SpeedProfile limits = build_speed_limit_profile_s( area, vehicle_params, comfort_settings, config, ego_v );
  if (!limits.empty()) {
      std::cout << "Speed limits built. s_max: " << limits.back().s << ", v_limit_max: " << limits.back().v_limit << std::endl;
  }

  // 3. ST Boundaries
  std::vector<StBoundary> boundaries;
  for( const auto& [id, p] : participants.participants )
  {
    auto b = build_boundary_for_participant( area, p, vehicle_params, config, ego_s, 0.0 );
    if( b )
    {
      boundaries.push_back( *b );
    }
  }

  // 4. Time samples
  const auto t_samples = build_time_samples( config.total_time, config.dt_dp );

  // 5. Build ST Forbidden Intervals
  auto forbidden = build_forbidden_intervals_per_t( t_samples, boundaries );
  std::cout << "ST Boundaries: " << boundaries.size() << ", forbidden intervals: " << forbidden.size() << std::endl;

  // 6. DP
  SpeedProfile dp_trajectory = solve_dp( config, comfort_settings, ego_s, ego_v, 0.0 /* ego_a? if avail */, 
                                         t_samples, forbidden, limits );
  
  if (dp_trajectory.empty()) {
      return {};
  }

  // 7. Tunnel & QP
  SpeedProfile final_trajectory = smooth_qp( config, ego_s, ego_v, 0.0, dp_trajectory, forbidden );

  // Fill in kappa, v_limit for output
  for( auto& p : final_trajectory )
  {
    p.kappa   = kappa_at_s_nearest( limits, p.s ); // approximations
    p.v_limit = v_limit_at_s_min_bracket( limits, p.s, config.default_legal_speed ); 
  }

  return final_trajectory;
}

} // namespace adore::planner
