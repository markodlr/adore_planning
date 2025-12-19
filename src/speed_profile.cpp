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

#include "adore_map/map_point.hpp"
#include "adore_math/curvature.hpp"
#include "adore_math/distance.h"
#include "adore_math/point.h"

#include <OsqpEigen/OsqpEigen.h>

namespace adore::planner
{
namespace
{

constexpr double k_inf = 1e30;

// =============================================================================
// Small helpers
// =============================================================================

static std::vector<double>
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

static double
interp_linear_monotone( const std::vector<double>& x, const std::vector<double>& y, double xq )
{
  const size_t n = x.size();
  if( n == 0 || y.size() != n )
  {
    return 0.0;
  }

  if( xq <= x.front() )
  {
    return y.front();
  }
  if( xq >= x.back() )
  {
    return y.back();
  }

  const auto   it_hi = std::lower_bound( x.begin(), x.end(), xq );
  const size_t i1    = static_cast<size_t>( it_hi - x.begin() );
  const size_t i0    = i1 - 1;

  const double x0 = x[i0];
  const double x1 = x[i1];
  const double d  = x1 - x0;
  const double u  = ( std::abs( d ) < 1e-12 ) ? 0.0 : ( xq - x0 ) / d;

  return y[i0] + std::clamp( u, 0.0, 1.0 ) * ( y[i1] - y[i0] );
}

static void
relax_monotone_bounds( std::vector<double>& lo, std::vector<double>& hi, double slack )
{
  for( size_t k = 1; k < lo.size(); ++k )
  {
    lo[k] = std::max( lo[k], lo[k - 1] - slack );
    hi[k] = std::max( hi[k], hi[k - 1] - slack );
  }
}

// =============================================================================
// Projection to reference line (ported)
// =============================================================================

struct SlPoint
{
  double s = 0.0;
  double l = 0.0;
};

static std::optional<SlPoint>
project_to_reference_line( const std::map<double, adore::map::MapPoint>& reference_line, const adore::math::Point2d& p,
                           std::optional<double> s_seed, double window )
{
  if( reference_line.size() < 2 )
  {
    return std::nullopt;
  }

  auto clamp01 = []( double u ) { return std::max( 0.0, std::min( 1.0, u ) ); };

  const double s_min = reference_line.begin()->first;
  const double s_max = reference_line.rbegin()->first;

  double win_lo = s_min;
  double win_hi = s_max;

  if( s_seed.has_value() && std::isfinite( *s_seed ) )
  {
    win_lo = std::max( s_min, *s_seed - window );
    win_hi = std::min( s_max, *s_seed + window );
    if( win_lo >= win_hi )
    {
      return std::nullopt;
    }
  }

  auto it0 = reference_line.lower_bound( win_lo );
  if( it0 == reference_line.end() )
  {
    return std::nullopt;
  }
  if( it0 != reference_line.begin() )
  {
    it0 = std::prev( it0 );
  }

  auto it_end = reference_line.upper_bound( win_hi );
  if( it_end == reference_line.begin() )
  {
    return std::nullopt;
  }

  double  best_d2 = std::numeric_limits<double>::infinity();
  SlPoint best;
  bool    ok = false;

  auto it      = it0;
  auto it_next = std::next( it );
  for( ; it_next != reference_line.end() && it_next != it_end; ++it, ++it_next )
  {
    const double s_a = it->first;
    const double s_b = it_next->first;

    const auto& a = it->second;
    const auto& b = it_next->second;

    const double vx       = b.x - a.x;
    const double vy       = b.y - a.y;
    const double seg_len2 = vx * vx + vy * vy;
    if( seg_len2 < 1e-9 )
    {
      continue;
    }

    const double wx = p.x - a.x;
    const double wy = p.y - a.y;

    const double u      = clamp01( ( wx * vx + wy * vy ) / seg_len2 );
    const double proj_x = a.x + u * vx;
    const double proj_y = a.y + u * vy;

    const double dx = p.x - proj_x;
    const double dy = p.y - proj_y;
    const double d2 = dx * dx + dy * dy;
    if( d2 >= best_d2 )
    {
      continue;
    }

    const double s_proj = s_a + u * ( s_b - s_a );

    const double seg_len = std::sqrt( seg_len2 );
    const double tx      = vx / seg_len;
    const double ty      = vy / seg_len;

    const double nx = -ty;
    const double ny = tx;

    const double l = ( p.x - proj_x ) * nx + ( p.y - proj_y ) * ny;

    best_d2 = d2;
    best.s  = s_proj;
    best.l  = l;
    ok      = true;
  }

  if( !ok )
  {
    return std::nullopt;
  }

  return best;
}

static std::optional<double>
ego_s_from_area( const DrivableArea& area, const dynamics::VehicleStateDynamic& ego, const SpeedProfileConfig& cfg )
{
  adore::math::Point2d p;
  p.x = ego.x;
  p.y = ego.y;

  const auto sl = project_to_reference_line( area.reference_line, p, std::nullopt, cfg.projection_window );
  if( !sl )
  {
    return std::nullopt;
  }
  return sl->s;
}

// =============================================================================
// Curvature + legal speed -> reference SpeedProfile along s
// =============================================================================

static void
smooth_kappa_in_place( std::vector<SpeedProfilePoint>& profile, size_t half_win )
{
  if( profile.size() < 3 || half_win <= 0 )
  {
    return;
  }
  const size_t n = profile.size();
  // Prefix sum of kappa (one extra element)
  std::vector<double> prefix( n + 1u, 0.0 );
  for( size_t i = 0; i < n; ++i )
  {
    prefix[i + 1u] = prefix[i] + profile[i].kappa;
  }

  for( size_t i = 0; i < n; ++i )
  {
    const size_t lo    = std::max( 0uz, i - half_win );
    const size_t hi    = std::min( n - 1, i + half_win );
    const size_t count = hi - lo + 1;
    const double sum   = prefix[hi + 1u] - prefix[lo];
    profile[i].kappa   = sum / static_cast<double>( count );
  }
}

static void
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

  // (Optional) recompute times using v (same as your existing time code, but use v not v_limit)
  profile.front().t = 0.0;
  for( size_t i = 0; i + 1 < profile.size(); ++i )
  {
    const double v_i   = std::max( 0.0, clamp_v( profile[i].v ) );
    const double v_ip1 = std::max( 0.0, clamp_v( profile[i + 1].v ) );
    const double dt    = ( std::abs( v_i + v_ip1 ) < 1e-6 ) ? 0.0 : ( 2.0 * ds ) / ( v_i + v_ip1 );
    profile[i + 1].t   = profile[i].t + dt;
  }
}

static SpeedProfile
build_speed_limit_profile_s( const DrivableArea& area, const SpeedProfileConfig& cfg, double v0 )
{
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
    if( kappa > cfg.max_curvature * 3.0 )
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
    v_legal        = std::max( 0.0, v_legal ) * cfg.speed_fraction_of_limit;
    v_legal        = std::clamp( v_legal, std::max( 0.0, cfg.v_min ), cfg.v_max );

    if( i == 0 )
    {
      v_legal = std::clamp( v0, std::max( 0.0, cfg.v_min ), cfg.v_max );
    }

    SpeedProfilePoint p;
    p.s       = s;
    p.t       = 0.0;
    p.kappa   = kappa;   // raw
    p.v_limit = v_legal; // legal for now
    ref_profile.push_back( p );
  }

  // smooth_kappa_in_place( ref_profile, cfg.curvature_smoothing_half_win );

  for( auto& p : ref_profile )
  {
    if( p.kappa < cfg.curvature_eps )
    {
      p.kappa = 0.0;
      continue;
    }
    const double v_curv = std::sqrt( std::max( 0.0, cfg.max_lateral_acc ) / p.kappa );
    const double v      = std::min( p.v_limit, v_curv );
    p.v_limit           = v;
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


  speed_forward_backward_pass( ref_profile, ds, cfg.a_max, cfg.a_min, cfg.v_min, cfg.v_max, v0 );

  return ref_profile;
}

// Conservative vmax lookup from speed_limit_profile_s (min of bracketing v_limit samples)
static double
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

static double
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
// ST boundaries (ported)
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

struct Interval
{
  double lo = 0.0;
  double hi = 0.0;
};

static bool
overlaps( const Interval& a, const Interval& b )
{
  return a.lo < b.hi && b.lo < a.hi;
}

static void
merge_intervals_in_place( std::vector<Interval>& intervals )
{
  if( intervals.empty() )
  {
    return;
  }

  std::sort( intervals.begin(), intervals.end(), []( const Interval& a, const Interval& b ) { return a.lo < b.lo; } );

  size_t out = 0;
  for( const auto& cur : intervals )
  {
    if( cur.hi <= cur.lo )
    {
      continue;
    }

    if( out == 0 || intervals[out - 1].hi < cur.lo )
    {
      intervals[out++] = cur;
    }
    else
    {
      intervals[out - 1].hi = std::max( intervals[out - 1].hi, cur.hi );
    }
  }

  intervals.resize( out );
}

static StBoundaryType
classify_boundary_type( const dynamics::TrafficParticipant& p )
{
  return p.trajectory.has_value() ? StBoundaryType::DynamicObstacle : StBoundaryType::StaticObstacle;
}

static std::optional<StBoundary>
build_boundary_for_participant( const DrivableArea& area, const dynamics::TrafficParticipant& p, const SpeedProfileConfig& cfg,
                                double ego_s_seed, double t0_abs )
{
  StBoundary boundary;
  boundary.obstacle_id = p.id;
  boundary.type        = classify_boundary_type( p );

  const size_t n = std::max( 1.0, std::ceil( cfg.total_time / cfg.obstacle_dt ) );
  boundary.samples.reserve( n + 1u );

  for( size_t i = 0; i <= n; ++i )
  {
    const double t_rel = static_cast<double>( i ) * cfg.obstacle_dt;
    const double t_abs = t0_abs + t_rel;

    const auto poly = p.get_corners_at_t( t_abs, cfg.obstacle_longitudinal_buffer, cfg.obstacle_lateral_buffer );
    if( poly.points.size() < 4 )
    {
      continue;
    }

    double s_min = std::numeric_limits<double>::infinity();
    double s_max = -std::numeric_limits<double>::infinity();

    for( const auto& c : poly.points )
    {
      adore::math::Point2d cp{ c.x, c.y };
      const auto           sl = project_to_reference_line( area.reference_line, cp, ego_s_seed, cfg.projection_window );
      if( !sl )
      {
        continue;
      }
      s_min = std::min( s_min, sl->s );
      s_max = std::max( s_max, sl->s );
    }

    if( !std::isfinite( s_min ) || !std::isfinite( s_max ) )
    {
      continue;
    }

    boundary.samples.push_back( { t_rel, s_min, s_max } );
  }

  if( boundary.samples.size() < 2 )
  {
    return std::nullopt;
  }

  std::sort( boundary.samples.begin(), boundary.samples.end(),
             []( const StBoundarySample& a, const StBoundarySample& b ) { return a.t < b.t; } );

  return boundary;
}

static std::vector<std::vector<Interval>>
build_forbidden_intervals_per_t( const std::vector<double>& t_samples, const std::vector<StBoundary>& boundaries )
{
  std::vector<std::vector<Interval>> forb( t_samples.size() );

  for( size_t k = 0; k < t_samples.size(); ++k )
  {
    const double          t = t_samples[k];
    std::vector<Interval> intervals;
    intervals.reserve( boundaries.size() );

    for( const auto& b : boundaries )
    {
      if( b.samples.size() < 2 )
      {
        continue;
      }

      if( t < b.samples.front().t || t > b.samples.back().t )
      {
        continue;
      }

      auto it = std::lower_bound( b.samples.begin(), b.samples.end(), t, []( const StBoundarySample& s, double tv ) { return s.t < tv; } );

      if( it == b.samples.begin() )
      {
        intervals.push_back( { it->s_lower, it->s_upper } );
        continue;
      }
      if( it == b.samples.end() )
      {
        intervals.push_back( { b.samples.back().s_lower, b.samples.back().s_upper } );
        continue;
      }

      const auto& s1 = *it;
      const auto& s0 = *std::prev( it );

      const double dt = s1.t - s0.t;
      const double u  = ( std::abs( dt ) < 1e-9 ) ? 0.0 : ( t - s0.t ) / dt;

      const double lo = s0.s_lower + u * ( s1.s_lower - s0.s_lower );
      const double hi = s0.s_upper + u * ( s1.s_upper - s0.s_upper );

      intervals.push_back( { std::min( lo, hi ), std::max( lo, hi ) } );
    }

    merge_intervals_in_place( intervals );
    forb[k] = std::move( intervals );
  }

  return forb;
}

static bool
is_forbidden( const std::vector<Interval>& forb, double s )
{
  for( const auto& in : forb )
  {
    if( s >= in.lo && s <= in.hi )
    {
      return true;
    }
  }
  return false;
}

// =============================================================================
// DP solver (ported, uses v_limit profile directly)
// =============================================================================

struct DpCell
{
  double cost = k_inf;
  int    prev = -1;
  double v    = 0.0;
  double a    = 0.0;
  bool   ok   = false;
};

static std::vector<double>
dp_solve( const SpeedProfileConfig& cfg, double dt, const SpeedProfile& speed_limit_profile_s, const std::vector<double>& t,
          const std::vector<std::vector<Interval>>& forb, double s0, double v0_in, double a0_in )
{
  if( t.size() < 2 )
  {
    std::cerr << "dp_solve: invalid time samples" << std::endl;
    return {};
  }
  const size_t n = t.size() - 1;

  const size_t m = std::max( 2.0, std::floor( cfg.s_horizon / cfg.ds_dp ) );

  auto s_at = [&]( size_t j ) { return s0 + static_cast<double>( j ) * cfg.ds_dp; };

  double v0 = v0_in;
  double a0 = a0_in;
  if( !std::isfinite( v0 ) || std::abs( v0 ) > 1e3 )
  {
    v0 = 0.0;
  }
  if( !std::isfinite( a0 ) || std::abs( a0 ) > 1e3 )
  {
    a0 = 0.0;
  }

  v0 = std::clamp( v0, cfg.v_min, cfg.v_max );
  a0 = std::clamp( a0, cfg.a_min, cfg.a_max );

  std::vector<std::vector<DpCell>> dp( n + 1, std::vector<DpCell>( m ) );

  dp[0][0].cost = 0.0;
  dp[0][0].prev = -1;
  dp[0][0].v    = v0;
  dp[0][0].a    = a0;
  dp[0][0].ok   = true;

  const double v_cruise = ( std::isfinite( cfg.v_cruise ) && cfg.v_cruise > 0.0 ) ? cfg.v_cruise : cfg.v_max;

  for( size_t k = 1; k <= n; ++k )
  {
    for( size_t j = 0; j < m; ++j )
    {
      const double s = s_at( j );
      if( is_forbidden( forb[k], s ) )
      {
        continue;
      }

      double best_cost = k_inf;
      int    best_prev = -1;
      double best_v    = 0.0;
      double best_a    = 0.0;

      const double s_max_step = cfg.v_max * dt;
      const double step_limit = std::ceil( s_max_step / cfg.ds_dp );
      const size_t jp_min     = ( j > step_limit ) ? static_cast<size_t>( j - step_limit ) : 0;

      for( size_t jp = jp_min; jp <= j; ++jp )
      {
        const auto& prev = dp[k - 1][jp];
        if( !prev.ok )
        {
          continue;
        }

        const double sp = s_at( jp );
        const double ds = s - sp;
        if( ds < 0.0 )
        {
          continue;
        }

        const double v = ds / dt;
        if( v < cfg.v_min - 1e-6 )
        {
          continue;
        }

        const double vmax_here = std::min( cfg.v_max, v_limit_at_s_min_bracket( speed_limit_profile_s, s, cfg.v_max ) );
        const double v_excess  = std::max( 0.0, v - vmax_here );

        const double a = ( v - prev.v ) / dt;
        if( a < cfg.a_min - 1e-6 || a > cfg.a_max + 1e-6 )
        {
          continue;
        }

        const double jrk = ( a - prev.a ) / dt;
        if( std::abs( jrk ) > cfg.j_max + 1e-6 )
        {
          continue;
        }

        const double v_ref      = std::min( vmax_here, v_cruise );
        const double c_progress = -cfg.w_dp_progress * ds;
        const double dv         = v - v_ref;
        const double c_speed    = cfg.w_dp_speed * ( dv * dv );
        const double c_accel    = cfg.w_dp_accel * ( a * a );
        const double c_jerk     = cfg.w_dp_jerk * ( jrk * jrk );
        const double c_vlimit   = cfg.w_dp_limit_violation * ( v_excess * v_excess );
        const double total_cost = prev.cost + c_progress + c_speed + c_accel + c_jerk + c_vlimit;

        if( total_cost < best_cost )
        {
          best_cost = total_cost;
          best_prev = static_cast<int>( jp );
          best_v    = v;
          best_a    = a;
        }
      }

      if( best_prev >= 0 && best_cost < k_inf / 2.0 )
      {
        auto& cell = dp[k][j];
        cell.ok    = true;
        cell.cost  = best_cost;
        cell.prev  = best_prev;
        cell.v     = best_v;
        cell.a     = best_a;
      }
    }
  }

  int    best_j    = -1;
  double best_cost = k_inf;
  for( size_t j = 0; j < m; ++j )
  {
    const auto& cell = dp[n][j];
    if( cell.ok && cell.cost < best_cost )
    {
      best_cost = cell.cost;
      best_j    = static_cast<int>( j );
    }
  }
  if( best_j < 0 )
  {
    std::cerr << "dp_solve: no feasible solution found" << std::endl;
    return {};
  }

  std::vector<double> s_dp( n + 1, s0 );
  int                 curr_j = best_j;
  for( size_t k = n; k <= n; --k )
  {
    s_dp[k] = s_at( static_cast<size_t>( curr_j ) );
    if( k == 0 )
      break;

    curr_j = dp[k][static_cast<size_t>( curr_j )].prev;
    if( curr_j < 0 )
    {
      std::cerr << "dp_solve: backtrack failed" << std::endl;
      return {};
    }
  }

  return s_dp;
}

// =============================================================================
// Tunnel builder (ported, with tunnel_expand support)
// =============================================================================

static void
build_tunnel( const SpeedProfileConfig& cfg, const std::vector<std::vector<Interval>>& forb, const std::vector<double>& s_dp,
              std::vector<double>& s_min, std::vector<double>& s_max )
{
  const size_t n = s_dp.size();
  s_min.assign( n, 0.0 );
  s_max.assign( n, 0.0 );

  const double back  = cfg.tunnel_back + cfg.tunnel_expand;
  const double front = cfg.tunnel_front + cfg.tunnel_expand;

  constexpr double eps = 0.05;

  for( size_t k = 0; k < n; ++k )
  {
    double lo = s_dp[k] - back;
    double hi = s_dp[k] + front;

    for( const auto& in : forb[k] )
    {
      Interval tunnel{ lo, hi };
      if( !overlaps( tunnel, in ) )
      {
        continue;
      }

      const double center = 0.5 * ( in.lo + in.hi );
      if( s_dp[k] <= center )
      {
        hi = std::min( hi, in.lo - eps );
      }
      else
      {
        lo = std::max( lo, in.hi + eps );
      }
    }

    if( lo > hi )
    {
      const double mid = 0.5 * ( lo + hi );
      lo               = mid;
      hi               = mid;
    }

    s_min[k] = lo;
    s_max[k] = hi;
  }

  // keep bounds non-decreasing to avoid infeasible monotonic constraints
  relax_monotone_bounds( s_min, s_max, 1e-3 );
}

// =============================================================================
// DP->QP upsampling (ported)
// =============================================================================

static void
upsample_track_and_tunnel( const std::vector<double>& t_coarse, const std::vector<double>& s_dp_coarse,
                           const std::vector<double>& s_min_coarse, const std::vector<double>& s_max_coarse,
                           const std::vector<double>& t_fine, std::vector<double>& s_dp_fine, std::vector<double>& s_min_fine,
                           std::vector<double>& s_max_fine )
{
  s_dp_fine.resize( t_fine.size() );
  s_min_fine.resize( t_fine.size() );
  s_max_fine.resize( t_fine.size() );

  for( size_t k = 0; k < t_fine.size(); ++k )
  {
    const double tf = t_fine[k];
    s_dp_fine[k]    = interp_linear_monotone( t_coarse, s_dp_coarse, tf );
    s_min_fine[k]   = interp_linear_monotone( t_coarse, s_min_coarse, tf );
    s_max_fine[k]   = interp_linear_monotone( t_coarse, s_max_coarse, tf );
    if( s_min_fine[k] > s_max_fine[k] )
    {
      std::swap( s_min_fine[k], s_max_fine[k] );
    }
  }

  relax_monotone_bounds( s_min_fine, s_max_fine, 1e-3 );
}

// =============================================================================
// QP smoother (ported; speed constraints use v_limit profile directly)
// =============================================================================

static void
clamp_kinematics_in_place( const SpeedProfileConfig& cfg, SpeedProfile& profile )
{
  for( auto& p : profile )
  {
    p.v = std::min( std::max( p.v, cfg.v_min ), cfg.v_max );
    p.a = std::min( std::max( p.a, cfg.a_min ), cfg.a_max );
    p.j = std::min( std::max( p.j, -cfg.j_max ), cfg.j_max );
  }
}

static SpeedProfile
qp_smooth( const SpeedProfileConfig& cfg, double dt, const SpeedProfile& speed_limit_profile_s, const std::vector<double>& t,
           const std::vector<double>& s_ref, const std::vector<double>& s_min, const std::vector<double>& s_max, double s0, double v0,
           double a0 )
{
  (void) a0;

  const bool k_debug_qp = false;

  const int nv = static_cast<int>( t.size() );
  if( nv < 2 )
  {
    return {};
  }
  if( static_cast<int>( s_min.size() ) != nv || static_cast<int>( s_max.size() ) != nv || static_cast<int>( s_ref.size() ) != nv )
  {
    std::cerr << "qp_smooth: tunnel/ref size mismatch\n";
    return {};
  }

  // Variables: v[0..nv-1], s[0..nv-1], e[0..nv-1] where e >= 0 is overspeed slack.
  const int idx_v = 0;
  const int idx_s = nv;
  const int idx_e = 2 * nv;
  const int nvar  = 3 * nv;

  auto v_i = [&]( int k ) { return idx_v + k; };
  auto s_i = [&]( int k ) { return idx_s + k; };
  auto e_i = [&]( int k ) { return idx_e + k; };

  // Constraints:
  // (1) v hard bounds: v[k] in [max(0,v_min), v_max]
  // (2) accel bounds
  // (3) jerk bounds
  // (4) s tunnel bounds
  // (5) dynamics: s[0]=s0, s[k+1]-s[k]-dt*v[k]=0
  // (6) v[0]=v0
  // (7) soft map speed limit: v[k] - e[k] <= vmax_map[k]
  // (8) slack bounds: e[k] >= 0 (no upper bound)

  const int n_v_bounds  = nv;
  const int n_acc       = nv - 1;
  const int n_jerk      = std::max( 0, nv - 2 );
  const int n_s_bounds  = nv;
  const int n_dyn       = nv; // s0 equality + (nv-1) dynamics
  const int n_v0_eq     = 1;
  const int n_vmap_soft = nv;
  const int n_e_bounds  = nv;

  const int nc = n_v_bounds + n_acc + n_jerk + n_s_bounds + n_dyn + n_v0_eq + n_vmap_soft + n_e_bounds;

  Eigen::SparseMatrix<double> A( nc, nvar );
  Eigen::VectorXd             lb( nc ), ub( nc );
  lb.setConstant( -k_inf );
  ub.setConstant( k_inf );

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve( static_cast<size_t>( 12 * nc ) );

  int row = 0;

  // (1) hard bounds on v: v[k] in [max(0,v_min), v_max]
  for( int k = 0; k < nv; ++k, ++row )
  {
    triplets.emplace_back( row, v_i( k ), 1.0 );
    lb( row ) = std::max( 0.0, cfg.v_min );
    ub( row ) = cfg.v_max;
  }

  // (2) accel bounds: v[k+1]-v[k] in [a_min*dt, a_max*dt]
  for( int k = 0; k < nv - 1; ++k, ++row )
  {
    triplets.emplace_back( row, v_i( k ), -1.0 );
    triplets.emplace_back( row, v_i( k + 1 ), 1.0 );
    lb( row ) = cfg.a_min * dt;
    ub( row ) = cfg.a_max * dt;
  }

  // (3) jerk bounds: v[k+2]-2v[k+1]+v[k] in [-j_max*dt^2, j_max*dt^2]
  for( int k = 0; k < nv - 2; ++k, ++row )
  {
    triplets.emplace_back( row, v_i( k ), 1.0 );
    triplets.emplace_back( row, v_i( k + 1 ), -2.0 );
    triplets.emplace_back( row, v_i( k + 2 ), 1.0 );
    lb( row ) = -cfg.j_max * dt * dt;
    ub( row ) = cfg.j_max * dt * dt;
  }

  // (4) tunnel bounds on s
  for( int k = 0; k < nv; ++k, ++row )
  {
    triplets.emplace_back( row, s_i( k ), 1.0 );
    lb( row ) = s_min[static_cast<size_t>( k )];
    ub( row ) = s_max[static_cast<size_t>( k )];
  }

  // (5) dynamics: s[0] == s0
  {
    triplets.emplace_back( row, s_i( 0 ), 1.0 );
    lb( row ) = s0;
    ub( row ) = s0;
    ++row;
  }

  // s[k+1] - s[k] - dt*v[k] == 0
  for( int k = 0; k < nv - 1; ++k, ++row )
  {
    triplets.emplace_back( row, s_i( k + 1 ), 1.0 );
    triplets.emplace_back( row, s_i( k ), -1.0 );
    triplets.emplace_back( row, v_i( k ), -dt );
    lb( row ) = 0.0;
    ub( row ) = 0.0;
  }

  // (6) v[0] == v0 (hard)
  {
    const double v0_clamped = std::clamp( std::isfinite( v0 ) ? v0 : 0.0, std::max( 0.0, cfg.v_min ), cfg.v_max );
    triplets.emplace_back( row, v_i( 0 ), 1.0 );
    lb( row ) = v0_clamped;
    ub( row ) = v0_clamped;
    ++row;
  }

  // (7) soft map speed limit: v[k] - e[k] <= vmax_map[k]
  // Implement as: (-inf) <= v[k] - e[k] <= vmax_map[k]
  for( int k = 0; k < nv; ++k, ++row )
  {
    triplets.emplace_back( row, v_i( k ), 1.0 );
    triplets.emplace_back( row, e_i( k ), -1.0 );

    const double s_eval   = std::clamp( s_ref[static_cast<size_t>( k )], s_min[static_cast<size_t>( k )], s_max[static_cast<size_t>( k )] );
    const double vmax_map = std::min( cfg.v_max, v_limit_at_s_min_bracket( speed_limit_profile_s, s_eval, cfg.v_max ) );

    lb( row ) = -k_inf;
    ub( row ) = vmax_map;
  }

  // (8) slack bounds: e[k] >= 0
  for( int k = 0; k < nv; ++k, ++row )
  {
    triplets.emplace_back( row, e_i( k ), 1.0 );
    lb( row ) = 0.0;
    ub( row ) = k_inf;
  }

  A.setFromTriplets( triplets.begin(), triplets.end() );

  // ---- Objective (only v gets track/acc/jerk; add strong penalty on e^2)
  Eigen::SparseMatrix<double> H( nvar, nvar );
  Eigen::VectorXd             g = Eigen::VectorXd::Zero( nvar );

  auto add_H = [&]( int i, int j, double v ) { H.coeffRef( i, j ) += v; };

  std::vector<double> v_ref( static_cast<size_t>( nv ), 0.0 );
  for( int k = 0; k < nv; ++k )
  {
    double v_est = 0.0;
    if( k == 0 )
    {
      v_est = ( s_ref[1] - s_ref[0] ) / dt;
    }
    else if( k == nv - 1 )
    {
      v_est = ( s_ref[static_cast<size_t>( nv - 1 )] - s_ref[static_cast<size_t>( nv - 2 )] ) / dt;
    }
    else
    {
      v_est = ( s_ref[static_cast<size_t>( k + 1 )] - s_ref[static_cast<size_t>( k - 1 )] ) / ( 2.0 * dt );
    }
    v_ref[static_cast<size_t>( k )] = std::clamp( v_est, std::max( 0.0, cfg.v_min ), cfg.v_max );
  }

  const double w_track = cfg.w_qp_track_dp;
  const double w_acc   = cfg.w_qp_accel;
  const double w_jerk  = cfg.w_qp_jerk;

  // Add this to your config; for now pick something big:
  const double w_over = 1e5; // overspeed penalty weight (tune)

  if( w_track > 0.0 )
  {
    for( int k = 0; k < nv; ++k )
    {
      const int i = v_i( k );
      add_H( i, i, 2.0 * w_track );
      g( i ) += -2.0 * w_track * v_ref[static_cast<size_t>( k )];
    }
  }

  if( w_acc > 0.0 && nv >= 2 )
  {
    for( int k = 0; k < nv - 1; ++k )
    {
      const int i0 = v_i( k );
      const int i1 = v_i( k + 1 );
      add_H( i0, i0, 2.0 * w_acc );
      add_H( i1, i1, 2.0 * w_acc );
      add_H( i0, i1, -2.0 * w_acc );
      add_H( i1, i0, -2.0 * w_acc );
    }
  }

  if( w_jerk > 0.0 && nv >= 3 )
  {
    for( int k = 0; k < nv - 2; ++k )
    {
      const int i0 = v_i( k );
      const int i1 = v_i( k + 1 );
      const int i2 = v_i( k + 2 );

      constexpr double c0 = 1.0;
      constexpr double c1 = -2.0;
      constexpr double c2 = 1.0;

      const int    idx[3] = { i0, i1, i2 };
      const double c[3]   = { c0, c1, c2 };

      for( int a = 0; a < 3; ++a )
      {
        for( int b = 0; b < 3; ++b )
        {
          add_H( idx[a], idx[b], 2.0 * w_jerk * c[a] * c[b] );
        }
      }
    }
  }

  // Overspeed slack penalty: sum w_over * e[k]^2
  if( w_over > 0.0 )
  {
    for( int k = 0; k < nv; ++k )
    {
      const int i = e_i( k );
      add_H( i, i, 2.0 * w_over );
    }
  }

  H.makeCompressed();

  OsqpEigen::Solver solver;
  solver.settings()->setVerbosity( false );
  solver.settings()->setWarmStart( true );
  solver.settings()->setMaxIteration( 10000 );
  solver.settings()->setAbsoluteTolerance( 1e-2 );

  solver.data()->setNumberOfVariables( nvar );
  solver.data()->setNumberOfConstraints( nc );

  if( !solver.data()->setHessianMatrix( H ) || !solver.data()->setGradient( g ) || !solver.data()->setLinearConstraintsMatrix( A )
      || !solver.data()->setLowerBound( lb ) || !solver.data()->setUpperBound( ub ) )
  {
    std::cerr << "qp_smooth: failed to set OSQP problem data\n";
    return {};
  }

  if( !solver.initSolver() )
  {
    std::cerr << "qp_smooth: initSolver failed\n";
    return {};
  }

  const auto err = solver.solveProblem();
  if( err != OsqpEigen::ErrorExitFlag::NoError )
  {
    std::cerr << "qp_smooth: solveProblem failed, err=" << static_cast<int>( err ) << "\n";
    return {};
  }

  // IMPORTANT: still check status; do not trust solution on infeasible.
  const auto status = solver.getStatus();
  if( status != OsqpEigen::Status::Solved )
  {
    std::cerr << "qp_smooth: OSQP status not solved: " << static_cast<int>( status ) << "\n";
    return {};
  }

  const Eigen::VectorXd sol = solver.getSolution();
  if( sol.size() != nvar )
  {
    std::cerr << "qp_smooth: unexpected solution size\n";
    return {};
  }

  if( k_debug_qp )
  {
    std::cerr << "qp_smooth: max overspeed slack=" << sol.segment( idx_e, nv ).maxCoeff() << "\n";
  }

  SpeedProfile out;
  out.resize( static_cast<size_t>( nv ) );

  for( int k = 0; k < nv; ++k )
  {
    const double v_raw = sol( v_i( k ) );
    const double s_raw = sol( s_i( k ) );

    if( !std::isfinite( v_raw ) || !std::isfinite( s_raw ) )
    {
      std::cerr << "qp_smooth: non-finite solution at k=" << k << "\n";
      return {};
    }

    out[static_cast<size_t>( k )].t = t[static_cast<size_t>( k )];
    out[static_cast<size_t>( k )].s = s_raw;
    out[static_cast<size_t>( k )].v = std::clamp( std::max( 0.0, v_raw ), std::max( 0.0, cfg.v_min ), cfg.v_max );
  }

  // accel from v
  for( int k = 0; k < nv; ++k )
  {
    double a = 0.0;
    if( nv >= 2 )
    {
      if( k == 0 )
        a = ( out[1].v - out[0].v ) / dt;
      else if( k == nv - 1 )
        a = ( out[static_cast<size_t>( nv - 1 )].v - out[static_cast<size_t>( nv - 2 )].v ) / dt;
      else
        a = ( out[static_cast<size_t>( k + 1 )].v - out[static_cast<size_t>( k - 1 )].v ) / ( 2.0 * dt );
    }
    out[static_cast<size_t>( k )].a = a;
  }

  // jerk from accel
  for( int k = 0; k < nv; ++k )
  {
    double j = 0.0;
    if( nv >= 3 )
    {
      if( k == 0 )
        j = ( out[1].a - out[0].a ) / dt;
      else if( k == nv - 1 )
        j = ( out[static_cast<size_t>( nv - 1 )].a - out[static_cast<size_t>( nv - 2 )].a ) / dt;
      else
        j = ( out[static_cast<size_t>( k + 1 )].a - out[static_cast<size_t>( k - 1 )].a ) / ( 2.0 * dt );
    }
    out[static_cast<size_t>( k )].j = j;
  }

  clamp_kinematics_in_place( cfg, out );
  return out;
}


} // namespace

// =============================================================================
// Public entrypoint
// =============================================================================

SpeedProfile
plan_speed_profile( const DrivableArea& area, const dynamics::TrafficParticipantSet& participants,
                    const dynamics::VehicleStateDynamic& ego_state, const SpeedProfileConfig& cfg )
{
  if( area.reference_line.size() < 2 )
  {
    return {};
  }

  const auto ego_s = ego_s_from_area( area, ego_state, cfg );
  if( !ego_s )
  {
    return {};
  }

  const double s0     = *ego_s;
  const double t0_abs = ego_state.time;

  // Speed limit profile along s (map legal + curvature -> v_limit)
  const SpeedProfile speed_limit_profile_s = build_speed_limit_profile_s( area, cfg, ego_state.vx );

  double min_v          = 13.6;
  double min_s          = std::numeric_limits<double>::quiet_NaN();
  double kappa_at_min_v = 0.0;
  for( const auto& p : speed_limit_profile_s )
  {
    if( p.v_limit < min_v )
    {
      min_v          = p.v_limit;
      min_s          = p.s;
      kappa_at_min_v = p.kappa;
    }
  }


  // time grids
  const auto t_dp = build_time_samples( cfg.total_time, cfg.dt_dp );
  const auto t_qp = build_time_samples( cfg.total_time, cfg.dt_qp );

  // boundaries
  std::vector<StBoundary> boundaries;
  boundaries.reserve( participants.participants.size() );

  for( const auto& [id, p] : participants.participants )
  {
    (void) id;
    auto b = build_boundary_for_participant( area, p, cfg, s0, t0_abs );
    if( b )
    {
      boundaries.push_back( std::move( *b ) );
    }
  }

  // forbidden intervals per DP time sample
  const auto forb_dp = build_forbidden_intervals_per_t( t_dp, boundaries );

  // DP solve
  const double v0 = ego_state.vx;
  const double a0 = ego_state.ax;


  const std::vector<double> s_dp = dp_solve( cfg, cfg.dt_dp, speed_limit_profile_s, t_dp, forb_dp, s0, v0, a0 );
  if( s_dp.empty() )
  {

    std::cerr << "v0=" << v0 << " a0=" << a0 << " v_min=" << cfg.v_min << " v_max=" << cfg.v_max << " a_min=" << cfg.a_min
              << " a_max=" << cfg.a_max << " j_max=" << cfg.j_max << "\n";

    const double v_t    = min_v;
    const double v0c    = v0;
    const double a      = std::abs( cfg.a_min );
    const double d_need = ( v0c * v0c - v_t * v_t ) / ( 2.0 * a );
    std::cerr << "brake_distance_needed=" << d_need << "  delta_to_min=" << ( min_s - s0 ) << "\n";

    std::cerr << "!!!!!!!!!!!!! Speed profile DP failed" << std::endl;
    return {};
  }

  // tunnel on DP grid
  std::vector<double> s_min_dp;
  std::vector<double> s_max_dp;
  build_tunnel( cfg, forb_dp, s_dp, s_min_dp, s_max_dp );

  // upsample DP + tunnel to QP grid
  std::vector<double> s_ref_qp;
  std::vector<double> s_min_qp;
  std::vector<double> s_max_qp;
  upsample_track_and_tunnel( t_dp, s_dp, s_min_dp, s_max_dp, t_qp, s_ref_qp, s_min_qp, s_max_qp );

  // QP smooth
  SpeedProfile out = qp_smooth( cfg, cfg.dt_qp, speed_limit_profile_s, t_qp, s_ref_qp, s_min_qp, s_max_qp, s0, v0, a0 );
  if( out.empty() )
  {
    std::cerr << "!!!!!!!!! Speed profile QP failed" << std::endl;
    return {};
  }

  // attach annotations
  for( auto& p : out )
  {
    p.v_limit = v_limit_at_s_min_bracket( speed_limit_profile_s, p.s, cfg.v_max );
    p.kappa   = kappa_at_s_nearest( speed_limit_profile_s, p.s );
  }

  // defensive clamp
  clamp_kinematics_in_place( cfg, out );
  return out;
}

} // namespace adore::planner
