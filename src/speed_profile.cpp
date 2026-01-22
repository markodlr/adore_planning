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

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "adore_map/map_point.hpp"
#include "adore_math/curvature.hpp"
#include "adore_math/distance.h"
#include "adore_math/geometry/interval.hpp"
#include "adore_math/geometry/projection.hpp"
#include "adore_math/point.h"

#include <OsqpEigen/OsqpEigen.h>

namespace adore::planner
{
namespace
{

constexpr double kInf                     = 1e30;
constexpr double kEpsilon                 = 1e-3;
constexpr double kBoundarySampleTolerance = 0.1;
constexpr double kStaticObstacleThreshold = 0.1; // [m/s]

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

StBoundaryType
classify_boundary_type( const dynamics::TrafficParticipant& p )
{
  const double v = std::hypot( p.state.vx, p.state.vy );
  // Simple heuristic: if slow, static.
  if( v < kStaticObstacleThreshold )
  {
    return StBoundaryType::StaticObstacle;
  }
  return StBoundaryType::DynamicObstacle;
}

std::optional<StBoundary>
build_boundary_for_participant( const DrivableArea& area, const dynamics::TrafficParticipant& p,
                                const dynamics::PhysicalVehicleParameters& vehicle_params, const SpeedProfileConfig& cfg, double ego_s_seed,
                                double t0_abs )
{
  StBoundary boundary;
  boundary.obstacle_id = p.id;
  boundary.type        = classify_boundary_type( p );

  const double t_horizon = cfg.total_time;
  const double dt        = cfg.obstacle_dt;
  const int    n_samples = static_cast<int>( std::ceil( t_horizon / dt ) );

  bool any_overlap = false;

  for( int i = 0; i <= n_samples; ++i )
  {
    const double t_rel = i * dt;
    const double v     = std::hypot( p.state.vx, p.state.vy );
    const double yaw   = p.state.yaw_angle;

    const double px = p.state.x + v * std::cos( yaw ) * t_rel;
    const double py = p.state.y + v * std::sin( yaw ) * t_rel;

    const auto& corners0 = p.get_corners().points;
    if( corners0.empty() )
      continue;

    double s_min = kInf, s_max = -kInf, l_min = kInf, l_max = -kInf;

    for( const auto& c0 : corners0 )
    {
      const double cx = px + ( c0.x - p.state.x );
      const double cy = py + ( c0.y - p.state.y );

      if( auto sl = adore::math::project_to_path( area.reference_line, adore::math::Point2d{ cx, cy }, std::optional<double>( ego_s_seed ),
                                                  cfg.projection_window ) )
      {
        s_min = std::min( s_min, sl->s );
        s_max = std::max( s_max, sl->s );
        l_min = std::min( l_min, sl->l );
        l_max = std::max( l_max, sl->l );
      }
    }

    if( s_min > s_max )
      continue;

    // Inflate
    s_min -= cfg.obstacle_longitudinal_buffer;
    s_max += cfg.obstacle_longitudinal_buffer;
    l_min -= cfg.obstacle_lateral_buffer;
    l_max += cfg.obstacle_lateral_buffer;

    const double s_avg    = 0.5 * ( s_min + s_max );
    auto         it_left  = area.left_boundary.lower_bound( s_avg );
    auto         it_right = area.right_boundary.lower_bound( s_avg );

    if( it_left == area.left_boundary.end() || it_right == area.right_boundary.end() )
      continue;

    auto sl_left  = adore::math::project_to_path( area.reference_line, adore::math::Point2d{ it_left->second.x, it_left->second.y },
                                                  std::optional<double>( s_avg ), 10.0 );
    auto sl_right = adore::math::project_to_path( area.reference_line, adore::math::Point2d{ it_right->second.x, it_right->second.y },
                                                  std::optional<double>( s_avg ), 10.0 );

    if( !sl_left || !sl_right )
      continue;

    const double safety_margin = ( vehicle_params.body_width * 0.5 ) + cfg.obstacle_lateral_buffer;
    const double overlap_min   = std::max( l_min, -safety_margin );
    const double overlap_max   = std::min( l_max, safety_margin );

    if( overlap_min < overlap_max )
    {
      boundary.samples.push_back( { t_rel, s_min, s_max } );
      any_overlap = true;
    }
  }

  if( any_overlap )
    return boundary;
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
      if( b.samples.empty() )
        continue;

      auto it = std::lower_bound( b.samples.begin(), b.samples.end(), t,
                                  []( const StBoundarySample& s, double val ) { return s.t < val; } );

      double s_lo = 0, s_hi = 0;
      bool   valid = false;

      if( it == b.samples.begin() )
      {
        if( std::abs( it->t - t ) < kBoundarySampleTolerance )
        {
          s_lo  = it->s_lower;
          s_hi  = it->s_upper;
          valid = true;
        }
      }
      else if( it == b.samples.end() )
      {
        auto prev = std::prev( it );
        if( std::abs( prev->t - t ) < kBoundarySampleTolerance )
        {
          s_lo  = prev->s_lower;
          s_hi  = prev->s_upper;
          valid = true;
        }
      }
      else
      {
        auto   prev = std::prev( it );
        double r    = ( t - prev->t ) / ( it->t - prev->t );
        s_lo        = prev->s_lower + r * ( it->s_lower - prev->s_lower );
        s_hi        = prev->s_upper + r * ( it->s_upper - prev->s_upper );
        valid       = true;
      }

      if( valid )
      {
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
  double cost = kInf;
  int    prev = -1;
  double v    = 0.0;
  double a    = 0.0;
  bool   ok   = false;
};

double
obstacle_proximity_cost( const SpeedProfileConfig& cfg, double s, double t, const std::vector<adore::math::Interval<double>>& forbidden )
{
  double min_dist = kInf;
  for( const auto& inv : forbidden )
  {
    double dist = 0.0;
    if( s < inv.lo )
      dist = inv.lo - s;
    else if( s > inv.hi )
      dist = s - inv.hi;
    else
      dist = 0.0;
    min_dist = std::min( min_dist, dist );
  }

  if( min_dist < kEpsilon )
    return cfg.w_dp_obstacle_proximity * 100.0;

  return cfg.w_dp_obstacle_proximity * std::exp( -min_dist );
}

SpeedProfile
solve_dp( const SpeedProfileConfig& cfg, const dynamics::ComfortSettings& comfort, double ego_s, double ego_v, double ego_a,
          const std::vector<double>& t_samples, const std::vector<std::vector<adore::math::Interval<double>>>& forbidden_intervals,
          const SpeedProfile& v_limit_profile, double s_ref_end )
{
  const size_t N_t = t_samples.size();

  if( N_t < 2 )
    return {};

  const double s_start = ego_s;
  // Clamp s_end to the end of the reference line
  const double s_end = std::min( ego_s + cfg.s_horizon, s_ref_end );
  const double ds    = cfg.ds_dp;
  const int    N_s   = static_cast<int>( std::ceil( ( s_end - s_start ) / ds ) );
  const int    width = N_s + 1;

  // Flattened cost table: [N_t * width]
  // Index: k * width + i
  std::vector<DpCell> cost_table( N_t * width );

  cost_table[0].cost = 0.0;
  cost_table[0].v    = ego_v;
  cost_table[0].a    = ego_a;
  cost_table[0].ok   = true;

  for( size_t k = 0; k < N_t - 1; ++k )
  {
    double dt = t_samples[k + 1] - t_samples[k];

    for( int i = 0; i <= N_s; ++i )
    {
      const int idx_curr = k * width + i;
      if( !cost_table[idx_curr].ok )
        continue;

      const double s_curr = s_start + i * ds;
      const double v_curr = cost_table[idx_curr].v;

      double s_next_min = s_curr + std::max( 0.0, v_curr * dt + 0.5 * comfort.min_acceleration * dt * dt );
      double s_next_max = s_curr + std::max( 0.0, v_curr * dt + 0.5 * comfort.max_acceleration * dt * dt );

      int j_min = static_cast<int>( ( s_next_min - s_start ) / ds );
      int j_max = static_cast<int>( ( s_next_max - s_start ) / ds ) + 1;

      j_min = std::max( i, j_min );
      j_max = std::min( N_s, j_max );

      for( int j = j_min; j <= j_max; ++j )
      {
        const int idx_next = ( k + 1 ) * width + j;

        double s_next = s_start + j * ds;
        double dist   = s_next - s_curr;
        double v_next = dist / dt * 2 - v_curr;
        double a      = ( v_next - v_curr ) / dt;

        const double local_v_limit = v_limit_at_s_min_bracket( v_limit_profile, s_next, comfort.max_speed );

        if( v_next < cfg.v_min - kEpsilon || v_next > comfort.max_speed + kEpsilon )
          continue;

        // Hard constraint: at or near the end of reference line, require v_next = 0
        const bool at_ref_end = ( s_next >= s_ref_end - ds - kEpsilon );
        if( at_ref_end && v_next > kEpsilon )
          continue;

        double v_excess = std::max( 0.0, v_next - local_v_limit );

        if( a < comfort.min_acceleration - kEpsilon || a > comfort.max_acceleration + kEpsilon )
          continue;

        bool collision = false;
        if( k + 1 < forbidden_intervals.size() )
        {
          for( const auto& inv : forbidden_intervals[k + 1] )
          {
            if( inv.contains( s_next ) )
            {
              collision = true;
              break;
            }
          }
        }
        if( collision )
          continue;

        double edge_cost = 0.0;
        double target_v  = ( cfg.v_cruise > 0.0 ) ? std::min( cfg.v_cruise, local_v_limit ) : local_v_limit;

        if( v_excess > kEpsilon )
        {
          edge_cost += cfg.w_dp_limit_violation * v_excess * v_excess;
        }

        edge_cost += cfg.w_dp_speed * std::abs( v_next - target_v );
        edge_cost += cfg.w_dp_accel * a * a;
        edge_cost += cfg.w_dp_jerk * std::abs( a - cost_table[idx_curr].a ) / dt;
        edge_cost += -cfg.w_dp_progress * dist;

        if( k + 1 < forbidden_intervals.size() )
          edge_cost += obstacle_proximity_cost( cfg, s_next, t_samples[k + 1], forbidden_intervals[k + 1] );

        double total = cost_table[idx_curr].cost + edge_cost;

        if( total < cost_table[idx_next].cost )
        {
          cost_table[idx_next].cost = total;
          cost_table[idx_next].prev = i;
          cost_table[idx_next].v    = v_next;
          cost_table[idx_next].a    = a;
          cost_table[idx_next].ok   = true;
        }
      }
    }
  }

  int       best_end_idx    = -1;
  double    min_c           = kInf;
  const int row_last_offset = ( N_t - 1 ) * width;

  for( int i = 0; i <= N_s; ++i )
  {
    int idx = row_last_offset + i;
    if( cost_table[idx].ok && cost_table[idx].cost < min_c )
    {
      min_c        = cost_table[idx].cost;
      best_end_idx = i;
    }
  }

  if( best_end_idx == -1 )
  {
    std::cerr << "solve_dp FAILED: no valid end state found" << std::endl;
    return {};
  }

  SpeedProfile result;
  result.resize( N_t );

  int curr_s_idx = best_end_idx;
  for( int k = N_t - 1; k >= 0; --k )
  {
    int idx = k * width + curr_s_idx;

    result[k].t = t_samples[k];
    result[k].s = s_start + curr_s_idx * ds;
    result[k].v = cost_table[idx].v;
    result[k].a = cost_table[idx].a;

    curr_s_idx = cost_table[idx].prev;
  }

  return result;
}

// =============================================================================
// Stop fence: forbid planning beyond reference line end
// =============================================================================

void
add_end_of_reference_fence( std::vector<std::vector<adore::math::Interval<double>>>& forbidden, double s_ref_end )
{
  // Add a forbidden interval [s_ref_end, +inf) for every time step
  for( auto& intervals_at_t : forbidden )
  {
    intervals_at_t.push_back( adore::math::Interval<double>( s_ref_end, kInf ) );
    adore::math::merge_intervals_in_place( intervals_at_t );
  }
}

// =============================================================================
// Fallback: brake-to-stop profile
// =============================================================================

SpeedProfile
build_brake_to_stop_profile( double ego_s, double ego_v, double s_ref_end, const dynamics::ComfortSettings& comfort,
                             const std::vector<double>& t_samples )
{
  SpeedProfile result;
  result.reserve( t_samples.size() );

  // Use comfort min_acceleration (negative) for braking
  const double a_brake = std::min( -0.1, comfort.min_acceleration ); // ensure negative

  double s = ego_s;
  double v = ego_v;

  for( size_t k = 0; k < t_samples.size(); ++k )
  {
    const double t  = t_samples[k];
    const double dt = ( k == 0 ) ? 0.0 : ( t_samples[k] - t_samples[k - 1] );

    if( k > 0 )
    {
      // Integrate kinematics
      const double v_new = std::max( 0.0, v + a_brake * dt );
      const double ds    = 0.5 * ( v + v_new ) * dt;
      s += ds;
      v = v_new;
    }

    // Clamp s to not exceed reference end
    s = std::min( s, s_ref_end );

    SpeedProfilePoint p;
    p.t       = t;
    p.s       = s;
    p.v       = v;
    p.a       = ( v > kEpsilon ) ? a_brake : 0.0;
    p.kappa   = 0.0;
    p.v_limit = 0.0;
    result.push_back( p );
  }

  return result;
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

  const auto sl = adore::math::project_to_path( area.reference_line, p, std::optional<double>( std::nullopt ), cfg.projection_window );
  if( !sl )
  {
    return std::nullopt;
  }
  return sl->s;
}

// =============================================================================
// Small helpers for main pipeline
// =============================================================================


SpeedProfile
build_speed_limit_profile_s( const DrivableArea& area, const dynamics::PhysicalVehicleParameters& vehicle_params,
                             const dynamics::ComfortSettings& comfort_settings, const SpeedProfileConfig& cfg, double v0 )
{
  // std::cerr << "build_speed_limit_profile_s: v0=" << v0 << ", comfort_v_max=" << comfort_settings.max_speed << std::endl;
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
    const double max_steer = std::max( std::abs( vehicle_params.steering_angle_max ), std::abs( vehicle_params.steering_angle_min ) );
    const double R_min     = vehicle_params.wheelbase / std::max( 0.01, std::tan( max_steer ) );
    const double physical_max_curvature  = 1.0 / std::max( 0.1, R_min );
    const double effective_max_curvature = std::min( cfg.max_curvature, physical_max_curvature );

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
    if( v_legal < 0.1 )
    {
      // std::cerr << " WARNING: v_legal=" << v_legal << " at s=" << s << std::endl;
    }
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
    // if( i % 10 == 0 )
    //   std::cerr << " i=" << i << " s=" << p.s << " v_legal=" << p.v_limit << " v_curv=" << v_curv << " final_v_limit=" << v << std::endl;
    p.v_limit = v;
    ++i;
  }
  const size_t n_tail = std::min<size_t>( 5u, ref_profile.size() );
  for( size_t k = 0; k < n_tail; ++k )
  {
    const size_t i         = ref_profile.size() - 1u - k;
    ref_profile[i].v_limit = 0.0;
  }


  return ref_profile;
}

} // namespace

// =============================================================================
// QP Smoother Logic (Soft Constraints)
// =============================================================================

namespace
{

// Resampled previous profile for QP temporal tracking
struct ResampledProfile
{
  std::vector<double> s;
  std::vector<double> v;
  std::vector<double> a;
  std::vector<bool>   valid;
};

// Build corridor bounds [s_lo, s_hi] for each time step
// Uses forbidden intervals + DP solution to define a convex tunnel
std::pair<std::vector<double>, std::vector<double>>
build_corridor_bounds_per_t( const std::vector<double>& t_samples, const std::vector<std::vector<adore::math::Interval<double>>>& forbidden,
                             const SpeedProfile& dp_profile, double ego_s, double s_ref_end, const SpeedProfileConfig& cfg )
{
  const size_t N = t_samples.size();

  std::vector<double> s_lo( N, ego_s );
  std::vector<double> s_hi( N, s_ref_end );

  for( size_t k = 0; k < N; ++k )
  {
    // Get DP solution at this time step
    double dp_s = ( k < dp_profile.size() ) ? dp_profile[k].s : s_ref_end;
    dp_s        = std::clamp( dp_s, ego_s, s_ref_end );

    // Compute free intervals as complement of forbidden within [ego_s, s_ref_end]
    std::vector<adore::math::Interval<double>> free_intervals;

    if( k < forbidden.size() && !forbidden[k].empty() )
    {
      // Sort forbidden intervals (they should already be merged)
      auto sorted_forbidden = forbidden[k];
      std::sort( sorted_forbidden.begin(), sorted_forbidden.end(), []( const auto& a, const auto& b ) { return a.lo < b.lo; } );

      double current = ego_s;
      for( const auto& fb : sorted_forbidden )
      {
        if( fb.lo > current )
        {
          // Free interval from current to fb.lo
          double lo = current + cfg.corridor_margin_s;
          double hi = fb.lo - cfg.corridor_margin_s;
          if( hi > lo && hi <= s_ref_end && lo >= ego_s )
          {
            free_intervals.push_back( { lo, hi } );
          }
        }
        current = std::max( current, fb.hi );
      }
      // Final free interval from last forbidden to s_ref_end
      if( current < s_ref_end )
      {
        double lo = current + cfg.corridor_margin_s;
        double hi = s_ref_end - cfg.corridor_margin_s;
        if( hi > lo )
        {
          free_intervals.push_back( { lo, hi } );
        }
      }
    }
    else
    {
      // No forbidden intervals - full range is free
      free_intervals.push_back( { ego_s, s_ref_end } );
    }

    // Select interval containing dp_s, or closest one
    bool   found     = false;
    double best_lo   = dp_s;
    double best_hi   = dp_s;
    double best_dist = kInf;

    for( const auto& inv : free_intervals )
    {
      if( inv.lo <= dp_s && dp_s <= inv.hi )
      {
        // DP point is inside this interval
        best_lo = inv.lo;
        best_hi = inv.hi;
        found   = true;
        break;
      }
      // Check distance to this interval
      double dist = 0.0;
      if( dp_s < inv.lo )
        dist = inv.lo - dp_s;
      else if( dp_s > inv.hi )
        dist = dp_s - inv.hi;

      if( dist < best_dist )
      {
        best_dist = dist;
        best_lo   = inv.lo;
        best_hi   = inv.hi;
      }
    }

    // Apply minimum width constraint
    if( best_hi - best_lo < cfg.corridor_min_width_s )
    {
      // Collapse around DP point
      double mid = std::clamp( dp_s, best_lo, best_hi );
      best_lo    = mid;
      best_hi    = mid;
    }

    // Enforce stop fence
    best_hi = std::min( best_hi, s_ref_end );
    best_lo = std::min( best_lo, s_ref_end );

    s_lo[k] = best_lo;
    s_hi[k] = best_hi;
  }

  return { s_lo, s_hi };
}

// Resample previous profile to current time grid for QP tracking
ResampledProfile
resample_previous_profile_to_t( const SpeedProfile& previous_profile, const std::vector<double>& t_samples )
{
  ResampledProfile result;
  const size_t     N = t_samples.size();
  result.s.resize( N, 0.0 );
  result.v.resize( N, 0.0 );
  result.a.resize( N, 0.0 );
  result.valid.resize( N, false );

  if( previous_profile.empty() )
  {
    return result;
  }

  for( size_t k = 0; k < N; ++k )
  {
    const double t_query = t_samples[k];

    // Find bracketing points in previous profile by time
    if( t_query < previous_profile.front().t || t_query > previous_profile.back().t )
    {
      continue; // Out of range
    }

    // Binary search for lower bound
    auto it = std::lower_bound( previous_profile.begin(), previous_profile.end(), t_query,
                                []( const SpeedProfilePoint& p, double t ) { return p.t < t; } );

    if( it == previous_profile.begin() )
    {
      result.s[k]     = it->s;
      result.v[k]     = it->v;
      result.a[k]     = it->a;
      result.valid[k] = true;
      continue;
    }

    auto prev_it = std::prev( it );
    if( it == previous_profile.end() )
    {
      result.s[k]     = prev_it->s;
      result.v[k]     = prev_it->v;
      result.a[k]     = prev_it->a;
      result.valid[k] = true;
      continue;
    }

    // Interpolate
    double t0    = prev_it->t;
    double t1    = it->t;
    double alpha = ( std::abs( t1 - t0 ) < 1e-9 ) ? 0.0 : ( t_query - t0 ) / ( t1 - t0 );
    alpha        = std::clamp( alpha, 0.0, 1.0 );

    result.s[k]     = prev_it->s + alpha * ( it->s - prev_it->s );
    result.v[k]     = prev_it->v + alpha * ( it->v - prev_it->v );
    result.a[k]     = prev_it->a + alpha * ( it->a - prev_it->a );
    result.valid[k] = true;
  }

  return result;
}

// Interpolate a value from a coarse time-indexed vector to a query time
double
interpolate_at_t( const std::vector<double>& t_coarse, const std::vector<double>& values, double t_query, double default_val )
{
  if( t_coarse.empty() || values.empty() )
    return default_val;

  if( t_query <= t_coarse.front() )
    return values.front();
  if( t_query >= t_coarse.back() )
    return values.back();

  // Binary search
  auto it = std::lower_bound( t_coarse.begin(), t_coarse.end(), t_query );
  if( it == t_coarse.begin() )
    return values.front();

  size_t idx_hi = static_cast<size_t>( std::distance( t_coarse.begin(), it ) );
  size_t idx_lo = idx_hi - 1;

  double t0    = t_coarse[idx_lo];
  double t1    = t_coarse[idx_hi];
  double alpha = ( std::abs( t1 - t0 ) < 1e-9 ) ? 0.0 : ( t_query - t0 ) / ( t1 - t0 );
  alpha        = std::clamp( alpha, 0.0, 1.0 );

  return values[idx_lo] + alpha * ( values[idx_hi] - values[idx_lo] );
}

// Resample DP profile to finer QP time grid
SpeedProfile
resample_dp_to_qp_grid( const SpeedProfile& dp_profile, const std::vector<double>& t_dp, const std::vector<double>& t_qp, double ego_s,
                        double ego_v, double ego_a )
{
  SpeedProfile result;
  result.reserve( t_qp.size() );

  if( dp_profile.empty() || t_dp.empty() )
    return result;

  // Build vectors for interpolation
  std::vector<double> dp_s( dp_profile.size() );
  std::vector<double> dp_v( dp_profile.size() );
  std::vector<double> dp_a( dp_profile.size() );
  for( size_t i = 0; i < dp_profile.size(); ++i )
  {
    dp_s[i] = dp_profile[i].s;
    dp_v[i] = dp_profile[i].v;
    dp_a[i] = dp_profile[i].a;
  }

  for( size_t k = 0; k < t_qp.size(); ++k )
  {
    double            t = t_qp[k];
    SpeedProfilePoint p;
    p.t = t;

    if( k == 0 )
    {
      // Force initial state to match ego exactly
      p.s = ego_s;
      p.v = ego_v;
      p.a = ego_a;
    }
    else
    {
      p.s = interpolate_at_t( t_dp, dp_s, t, ego_s );
      p.v = interpolate_at_t( t_dp, dp_v, t, 0.0 );
      p.a = interpolate_at_t( t_dp, dp_a, t, 0.0 );
    }

    result.push_back( p );
  }

  return result;
}

// Resample corridor bounds to finer QP time grid
std::pair<std::vector<double>, std::vector<double>>
resample_corridor_to_qp_grid( const std::vector<double>& s_lo_dp, const std::vector<double>& s_hi_dp, const std::vector<double>& t_dp,
                              const std::vector<double>& t_qp, double ego_s, double s_ref_end )
{
  const size_t        N_qp = t_qp.size();
  std::vector<double> s_lo_qp( N_qp );
  std::vector<double> s_hi_qp( N_qp );

  for( size_t k = 0; k < N_qp; ++k )
  {
    double t = t_qp[k];

    if( k == 0 )
    {
      // At t=0, corridor should allow starting from ego_s
      s_lo_qp[k] = ego_s;
      s_hi_qp[k] = s_ref_end;
    }
    else
    {
      s_lo_qp[k] = interpolate_at_t( t_dp, s_lo_dp, t, ego_s );
      s_hi_qp[k] = interpolate_at_t( t_dp, s_hi_dp, t, s_ref_end );
    }

    // Ensure bounds are valid
    s_lo_qp[k] = std::max( ego_s, s_lo_qp[k] );
    s_hi_qp[k] = std::min( s_ref_end, s_hi_qp[k] );
  }

  return { s_lo_qp, s_hi_qp };
}

// Main OSQP-based QP smoother
SpeedProfile
smooth_speed_profile_qp_osqp( const std::vector<double>& t_samples, const SpeedProfile& dp_profile, const std::vector<double>& s_lo,
                              const std::vector<double>& s_hi, const ResampledProfile& prev_resampled, double ego_s, double ego_v,
                              double ego_a, double s_ref_end, const dynamics::ComfortSettings& comfort, const SpeedProfileConfig& cfg )
{
  const size_t N = t_samples.size();
  if( N < 2 || dp_profile.size() < 2 )
  {
    return dp_profile; // Fallback
  }

  // Decision variables: [s(0..N-1), v(0..N-1), a(0..N-1), slack_lo(0..N-1), slack_hi(0..N-1)]
  const int n_vars = 5 * static_cast<int>( N );

  // Constraint count:
  // - N initial/dynamics constraints (s evolution): s[k+1] = s[k] + v[k]*dt + 0.5*a[k]*dt^2
  // - N-1 velocity evolution: v[k+1] = v[k] + a[k]*dt
  // - Initial conditions: s[0], v[0], a[0] = ego values (3 constraints)
  // - Corridor constraints with slacks: 2*N
  // - Slack nonnegativity: 2*N
  // - Monotonic station: N-1
  // - Speed bounds: N
  // - Accel bounds: N

  // For simplicity, encode most as box constraints on variables and use linear constraints for dynamics

  // Indices in decision vector:
  auto idx_s  = []( size_t k ) { return static_cast<int>( k ); };
  auto idx_v  = [N]( size_t k ) { return static_cast<int>( N + k ); };
  auto idx_a  = [N]( size_t k ) { return static_cast<int>( 2 * N + k ); };
  auto idx_sl = [N]( size_t k ) { return static_cast<int>( 3 * N + k ); }; // slack_lo
  auto idx_sh = [N]( size_t k ) { return static_cast<int>( 4 * N + k ); }; // slack_hi

  // Build Hessian (sparse, diagonal-dominant for this formulation)
  std::vector<Eigen::Triplet<double>> H_triplets;
  Eigen::VectorXd                     gradient = Eigen::VectorXd::Zero( n_vars );

  // Cost: track DP, track previous, accel, jerk, slack, terminal
  for( size_t k = 0; k < N; ++k )
  {
    // Track DP station
    double dp_s = ( k < dp_profile.size() ) ? dp_profile[k].s : s_ref_end;
    H_triplets.emplace_back( idx_s( k ), idx_s( k ), cfg.w_qp_track_dp );
    gradient( idx_s( k ) ) -= cfg.w_qp_track_dp * dp_s;

    // Track previous (if valid)
    if( k < prev_resampled.valid.size() && prev_resampled.valid[k] )
    {
      H_triplets.emplace_back( idx_s( k ), idx_s( k ), cfg.w_qp_track_prev );
      gradient( idx_s( k ) ) -= cfg.w_qp_track_prev * prev_resampled.s[k];
    }

    // Acceleration penalty
    H_triplets.emplace_back( idx_a( k ), idx_a( k ), cfg.w_qp_accel );

    // Slack penalties
    H_triplets.emplace_back( idx_sl( k ), idx_sl( k ), cfg.w_qp_slack );
    H_triplets.emplace_back( idx_sh( k ), idx_sh( k ), cfg.w_qp_slack );
  }

  // Jerk penalty: (a[k+1] - a[k])^2
  for( size_t k = 0; k + 1 < N; ++k )
  {
    // Expands to: a[k+1]^2 - 2*a[k+1]*a[k] + a[k]^2
    H_triplets.emplace_back( idx_a( k ), idx_a( k ), cfg.w_qp_jerk );
    H_triplets.emplace_back( idx_a( k + 1 ), idx_a( k + 1 ), cfg.w_qp_jerk );
    H_triplets.emplace_back( idx_a( k ), idx_a( k + 1 ), -cfg.w_qp_jerk );
    H_triplets.emplace_back( idx_a( k + 1 ), idx_a( k ), -cfg.w_qp_jerk );
  }

  // Terminal cost: penalize v[N-1]^2 and (s[N-1] - s_ref_end)^2
  H_triplets.emplace_back( idx_v( N - 1 ), idx_v( N - 1 ), cfg.w_qp_terminal );
  H_triplets.emplace_back( idx_s( N - 1 ), idx_s( N - 1 ), cfg.w_qp_terminal );
  gradient( idx_s( N - 1 ) ) -= cfg.w_qp_terminal * s_ref_end;

  // Build Hessian matrix
  Eigen::SparseMatrix<double> H( n_vars, n_vars );
  H.setFromTriplets( H_triplets.begin(), H_triplets.end() );
  H = ( H + Eigen::SparseMatrix<double>( H.transpose() ) ) / 2.0; // Ensure symmetry

  // Constraints:
  // Dynamics: s[k+1] - s[k] - v[k]*dt - 0.5*a[k]*dt^2 = 0  for k=0..N-2
  // Velocity: v[k+1] - v[k] - a[k]*dt = 0  for k=0..N-2
  // Initial: s[0] = ego_s, v[0] = ego_v, a[0] = ego_a
  // Corridor: s[k] + slack_lo[k] >= s_lo[k], s[k] - slack_hi[k] <= s_hi[k]
  // Monotonic: s[k+1] - s[k] >= 0

  const int n_dyn_s   = static_cast<int>( N - 1 );
  const int n_dyn_v   = static_cast<int>( N - 1 );
  const int n_init    = 3;
  const int n_corr_lo = static_cast<int>( N );
  const int n_corr_hi = static_cast<int>( N );
  const int n_mono    = static_cast<int>( N - 1 );
  const int n_cons    = n_dyn_s + n_dyn_v + n_init + n_corr_lo + n_corr_hi + n_mono;

  std::vector<Eigen::Triplet<double>> A_triplets;
  Eigen::VectorXd                     lb = Eigen::VectorXd::Constant( n_cons, -OsqpEigen::INFTY );
  Eigen::VectorXd                     ub = Eigen::VectorXd::Constant( n_cons, OsqpEigen::INFTY );

  int row = 0;

  // Dynamics constraints for s: s[k+1] - s[k] - v[k]*dt - 0.5*a[k]*dt^2 = 0
  for( size_t k = 0; k + 1 < N; ++k )
  {
    double dt = t_samples[k + 1] - t_samples[k];
    A_triplets.emplace_back( row, idx_s( k + 1 ), 1.0 );
    A_triplets.emplace_back( row, idx_s( k ), -1.0 );
    A_triplets.emplace_back( row, idx_v( k ), -dt );
    A_triplets.emplace_back( row, idx_a( k ), -0.5 * dt * dt );
    lb( row ) = 0.0;
    ub( row ) = 0.0;
    ++row;
  }

  // Dynamics constraints for v: v[k+1] - v[k] - a[k]*dt = 0
  for( size_t k = 0; k + 1 < N; ++k )
  {
    double dt = t_samples[k + 1] - t_samples[k];
    A_triplets.emplace_back( row, idx_v( k + 1 ), 1.0 );
    A_triplets.emplace_back( row, idx_v( k ), -1.0 );
    A_triplets.emplace_back( row, idx_a( k ), -dt );
    lb( row ) = 0.0;
    ub( row ) = 0.0;
    ++row;
  }

  // Initial conditions
  A_triplets.emplace_back( row, idx_s( 0 ), 1.0 );
  lb( row ) = ego_s;
  ub( row ) = ego_s;
  ++row;

  A_triplets.emplace_back( row, idx_v( 0 ), 1.0 );
  lb( row ) = ego_v;
  ub( row ) = ego_v;
  ++row;

  A_triplets.emplace_back( row, idx_a( 0 ), 1.0 );
  lb( row ) = ego_a;
  ub( row ) = ego_a;
  ++row;

  // Corridor constraints: s[k] + slack_lo[k] >= s_lo[k]
  for( size_t k = 0; k < N; ++k )
  {
    A_triplets.emplace_back( row, idx_s( k ), 1.0 );
    A_triplets.emplace_back( row, idx_sl( k ), 1.0 );
    lb( row ) = s_lo[k];
    ub( row ) = OsqpEigen::INFTY;
    ++row;
  }

  // Corridor constraints: s[k] - slack_hi[k] <= s_hi[k]
  for( size_t k = 0; k < N; ++k )
  {
    A_triplets.emplace_back( row, idx_s( k ), 1.0 );
    A_triplets.emplace_back( row, idx_sh( k ), -1.0 );
    lb( row ) = -OsqpEigen::INFTY;
    ub( row ) = s_hi[k];
    ++row;
  }

  // Monotonic constraints: s[k+1] - s[k] >= 0
  for( size_t k = 0; k + 1 < N; ++k )
  {
    A_triplets.emplace_back( row, idx_s( k + 1 ), 1.0 );
    A_triplets.emplace_back( row, idx_s( k ), -1.0 );
    lb( row ) = 0.0;
    ub( row ) = OsqpEigen::INFTY;
    ++row;
  }

  Eigen::SparseMatrix<double> A( n_cons, n_vars );
  A.setFromTriplets( A_triplets.begin(), A_triplets.end() );

  // Variable bounds (box constraints)
  Eigen::VectorXd x_lb = Eigen::VectorXd::Constant( n_vars, -OsqpEigen::INFTY );
  Eigen::VectorXd x_ub = Eigen::VectorXd::Constant( n_vars, OsqpEigen::INFTY );

  for( size_t k = 0; k < N; ++k )
  {
    // Station bounds (loose)
    x_lb( idx_s( k ) ) = ego_s;
    x_ub( idx_s( k ) ) = s_ref_end + cfg.max_slack_s;

    // Velocity bounds
    x_lb( idx_v( k ) ) = 0.0;
    x_ub( idx_v( k ) ) = comfort.max_speed;

    // Acceleration bounds
    x_lb( idx_a( k ) ) = comfort.min_acceleration;
    x_ub( idx_a( k ) ) = comfort.max_acceleration;

    // Slack bounds (nonnegative, capped)
    x_lb( idx_sl( k ) ) = 0.0;
    x_ub( idx_sl( k ) ) = cfg.max_slack_s;
    x_lb( idx_sh( k ) ) = 0.0;
    x_ub( idx_sh( k ) ) = cfg.max_slack_s;
  }

  // Combine linear constraints with box constraints
  // OSQP requires: l <= Ax <= u, so we add identity rows for box constraints
  const int                   n_box      = n_vars;
  const int                   total_cons = n_cons + n_box;
  Eigen::SparseMatrix<double> A_full( total_cons, n_vars );
  Eigen::VectorXd             lb_full( total_cons );
  Eigen::VectorXd             ub_full( total_cons );

  // Copy linear constraints
  std::vector<Eigen::Triplet<double>> A_full_triplets = A_triplets;
  lb_full.head( n_cons )                              = lb;
  ub_full.head( n_cons )                              = ub;

  // Add identity for box constraints
  for( int i = 0; i < n_vars; ++i )
  {
    A_full_triplets.emplace_back( n_cons + i, i, 1.0 );
    lb_full( n_cons + i ) = x_lb( i );
    ub_full( n_cons + i ) = x_ub( i );
  }
  A_full.setFromTriplets( A_full_triplets.begin(), A_full_triplets.end() );

  // Setup OSQP solver
  OsqpEigen::Solver solver;
  solver.settings()->setWarmStart( cfg.osqp_warm_start );
  solver.settings()->setVerbosity( false );
  solver.settings()->setPolish( cfg.osqp_polish );
  solver.settings()->setMaxIteration( cfg.osqp_max_iter );
  solver.settings()->setAbsoluteTolerance( cfg.osqp_eps_abs );
  solver.settings()->setRelativeTolerance( cfg.osqp_eps_rel );

  solver.data()->setNumberOfVariables( n_vars );
  solver.data()->setNumberOfConstraints( total_cons );

  if( !solver.data()->setHessianMatrix( H ) || !solver.data()->setGradient( gradient )
      || !solver.data()->setLinearConstraintsMatrix( A_full ) || !solver.data()->setLowerBound( lb_full )
      || !solver.data()->setUpperBound( ub_full ) )
  {
    return dp_profile; // Fallback
  }

  if( !solver.initSolver() )
  {
    return dp_profile; // Fallback
  }

  // Warm start from DP solution (primal variables only)
  Eigen::VectorXd x0 = Eigen::VectorXd::Zero( n_vars );
  for( size_t k = 0; k < N; ++k )
  {
    x0( idx_s( k ) )  = ( k < dp_profile.size() ) ? dp_profile[k].s : s_ref_end;
    x0( idx_v( k ) )  = ( k < dp_profile.size() ) ? dp_profile[k].v : 0.0;
    x0( idx_a( k ) )  = ( k < dp_profile.size() ) ? dp_profile[k].a : 0.0;
    x0( idx_sl( k ) ) = 0.0;
    x0( idx_sh( k ) ) = 0.0;
  }
  solver.setPrimalVariable( x0 );

  // Solve
  auto status = solver.solveProblem();
  if( status != OsqpEigen::ErrorExitFlag::NoError )
  {
    return dp_profile; // Fallback
  }

  // Extract solution
  Eigen::VectorXd solution = solver.getSolution();

  SpeedProfile result;
  result.reserve( N );

  for( size_t k = 0; k < N; ++k )
  {
    SpeedProfilePoint p;
    p.t = t_samples[k];
    p.s = std::clamp( solution( idx_s( k ) ), ego_s, s_ref_end );
    p.v = std::max( 0.0, solution( idx_v( k ) ) );
    p.a = solution( idx_a( k ) );
    result.push_back( p );
  }

  return result;
}

} // namespace

// =============================================================================
// Public API Implementation
// =============================================================================

SpeedProfile
plan_speed_profile( const DrivableArea& area, const dynamics::TrafficParticipantSet& participants,
                    const dynamics::VehicleStateDynamic& ego_state, const dynamics::PhysicalVehicleParameters& vehicle_params,
                    const dynamics::ComfortSettings& comfort_settings, const SpeedProfileConfig& config,
                    const SpeedProfile* previous_profile )
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
  // Use actual ego acceleration (clamped to comfort bounds for robustness)
  double ego_a = std::clamp( ego_state.ax, comfort_settings.min_acceleration, comfort_settings.max_acceleration );

  // 2. Speed Limits
  SpeedProfile limits = build_speed_limit_profile_s( area, vehicle_params, comfort_settings, config, ego_v );

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

  // 5b. Add stop fence: forbid s >= s_ref_end for all time steps
  const double s_ref_end = area.reference_line.rbegin()->first;
  add_end_of_reference_fence( forbidden, s_ref_end );

  // 6. DP - pass the reference line end for hard stop constraint
  SpeedProfile dp_trajectory = solve_dp( config, comfort_settings, ego_s, ego_v, ego_a, t_samples, forbidden, limits, s_ref_end );

  // 7. Fallback: if DP fails, use brake-to-stop profile
  if( dp_trajectory.empty() )
  {
    dp_trajectory = build_brake_to_stop_profile( ego_s, ego_v, s_ref_end, comfort_settings, t_samples );
  }

  SpeedProfile final_trajectory = dp_trajectory;

  // 8. QP smoothing (if enabled) - uses finer dt_qp resolution
  if( false && config.enable_qp_smoothing && !dp_trajectory.empty() )
  {
    // Build corridor bounds on coarse DP grid
    auto [s_lo_dp, s_hi_dp] = build_corridor_bounds_per_t( t_samples, forbidden, dp_trajectory, ego_s, s_ref_end, config );

    // Create finer QP time grid using dt_qp
    const auto t_qp = build_time_samples( config.total_time, config.dt_qp );

    // Resample DP profile to finer QP grid (forces ego state at t=0)
    SpeedProfile dp_resampled = resample_dp_to_qp_grid( dp_trajectory, t_samples, t_qp, ego_s, ego_v, ego_a );

    // Resample corridor bounds to finer QP grid
    auto [s_lo_qp, s_hi_qp] = resample_corridor_to_qp_grid( s_lo_dp, s_hi_dp, t_samples, t_qp, ego_s, s_ref_end );

    // Resample previous profile for temporal tracking (on QP grid)
    ResampledProfile prev_resampled;
    if( previous_profile && !previous_profile->empty() )
    {
      prev_resampled = resample_previous_profile_to_t( *previous_profile, t_qp );
    }

    // Run QP smoother on finer grid
    SpeedProfile qp_result = smooth_speed_profile_qp_osqp( t_qp, dp_resampled, s_lo_qp, s_hi_qp, prev_resampled, ego_s, ego_v, ego_a,
                                                           s_ref_end, comfort_settings, config );

    if( !qp_result.empty() )
    {
      final_trajectory = qp_result;
    }
  }

  // Fill in kappa, v_limit for output
  for( auto& p : final_trajectory )
  {
    p.kappa   = kappa_at_s_nearest( limits, p.s );
    p.v_limit = v_limit_at_s_min_bracket( limits, p.s, config.default_legal_speed );
  }

  return final_trajectory;
}


} // namespace adore::planner
