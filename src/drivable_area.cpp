/********************************************************************************
 * Copyright (c) 2025 Contributors
 *
 * SPDX-License-Identifier: EPL-2.0
 ********************************************************************************/

#include "planning/drivable_area.hpp"

#include <cmath>

#include <chrono>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <OsqpEigen/OsqpEigen.h>

namespace adore
{
namespace planner
{

namespace
{
// =============================================================================
// constants
// =============================================================================
constexpr double k_min_width           = 1.0;  // [m]
constexpr double k_station_epsilon     = 1e-3; // [m]
constexpr double k_small_long_buffer   = 0.5;  // [m]
constexpr double k_heading_probe_ds    = 0.1;  // [m]
constexpr double k_participant_s_range = 10.0;

// =============================================================================
// small math helpers
// =============================================================================
inline double
dot2( double ax, double ay, double bx, double by )
{
  return ax * bx + ay * by;
}

inline double
clamp01( double u )
{
  return std::max( 0.0, std::min( 1.0, u ) );
}

struct SlPoint
{
  double s  = 0.0; // drivable-area station (key of reference_line)
  double l  = 0.0; // signed lateral offset (+left of travel)
  bool   ok = false;
};

inline std::optional<SlPoint>
project_to_reference_line( const std::map<double, adore::map::MapPoint>& ref_line, const adore::math::Point2d& p,
                           std::optional<double> s_seed, double max_search )
{
  if( ref_line.size() < 2 )
    return std::nullopt;

  const double s_min = ref_line.begin()->first;
  const double s_max = ref_line.rbegin()->first;

  double win_lo = s_min;
  double win_hi = s_max;

  if( s_seed.has_value() && std::isfinite( *s_seed ) )
  {
    win_lo = std::max( s_min, *s_seed - max_search );
    win_hi = std::min( s_max, *s_seed + max_search );
    if( win_lo >= win_hi )
      return std::nullopt;
  }

  auto it0 = ref_line.lower_bound( win_lo );
  if( it0 == ref_line.end() )
    return std::nullopt;
  if( it0 != ref_line.begin() )
    it0 = std::prev( it0 ); // include segment crossing window start

  auto it_end = ref_line.upper_bound( win_hi );
  if( it_end == ref_line.begin() )
    return std::nullopt;

  double  best_d2 = std::numeric_limits<double>::infinity();
  SlPoint best;

  auto it      = it0;
  auto it_next = std::next( it );
  for( ; it_next != ref_line.end() && it_next != it_end; ++it, ++it_next )
  {
    const double s_a = it->first;
    const double s_b = it_next->first;

    const auto& a = it->second;
    const auto& b = it_next->second;

    const double vx       = b.x - a.x;
    const double vy       = b.y - a.y;
    const double seg_len2 = vx * vx + vy * vy;
    if( seg_len2 < 1e-9 )
      continue;

    const double wx = p.x - a.x;
    const double wy = p.y - a.y;

    const double u      = clamp01( ( wx * vx + wy * vy ) / seg_len2 );
    const double proj_x = a.x + u * vx;
    const double proj_y = a.y + u * vy;

    const double dx = p.x - proj_x;
    const double dy = p.y - proj_y;
    const double d2 = dx * dx + dy * dy;

    if( d2 >= best_d2 )
      continue;

    const double s_proj = s_a + u * ( s_b - s_a );

    const double seg_len = std::sqrt( seg_len2 );
    const double tx      = vx / seg_len;
    const double ty      = vy / seg_len;

    // left normal (Frenet +l)
    const double nx = -ty;
    const double ny = tx;

    const double l = ( p.x - proj_x ) * nx + ( p.y - proj_y ) * ny;

    best_d2 = d2;
    best.s  = s_proj;
    best.l  = l;
    best.ok = true;
  }

  if( !best.ok )
    return std::nullopt;

  return best;
}

// =============================================================================
// Route-aligned local frame
// =============================================================================
struct RouteFrame
{
  double cx = 0.0;
  double cy = 0.0;
  double nx = 0.0; // left normal
  double ny = 0.0;
  bool   ok = false;
};

inline RouteFrame
make_route_frame( const adore::map::Route& route, double s )
{
  const auto p  = route.get_map_point_at_s( s );
  const auto p2 = route.get_map_point_at_s( s + k_heading_probe_ds );

  const double dx  = p2.x - p.x;
  const double dy  = p2.y - p.y;
  const double len = std::hypot( dx, dy );

  const double angle = std::atan2( dy, dx );
  return { p.x, p.y, -std::sin( angle ), std::cos( angle ), true };
}

inline double
lateral_offset( const RouteFrame& f, double x, double y )
{
  return ( x - f.cx ) * f.nx + ( y - f.cy ) * f.ny;
}

// =============================================================================
// Oriented rectangle (for participant carving)
// =============================================================================
struct OrientedRect
{
  double cx = 0.0;
  double cy = 0.0;

  double ux = 1.0;
  double uy = 0.0;
  double vx = 0.0;
  double vy = 1.0;

  double hu = 0.0;
  double hv = 0.0;
};

bool
make_oriented_rect_from_corners( const std::vector<math::Point2d>& corners, OrientedRect& r )
{
  if( corners.size() < 4 )
    return false;

  // center
  r.cx = r.cy = 0.0;
  for( int i = 0; i < 4; ++i )
  {
    r.cx += corners[i].x;
    r.cy += corners[i].y;
  }
  r.cx *= 0.25;
  r.cy *= 0.25;

  // covariance
  double sxx = 0.0, sxy = 0.0, syy = 0.0;
  for( int i = 0; i < 4; ++i )
  {
    const double dx  = corners[i].x - r.cx;
    const double dy  = corners[i].y - r.cy;
    sxx             += dx * dx;
    sxy             += dx * dy;
    syy             += dy * dy;
  }
  sxx *= 0.25;
  sxy *= 0.25;
  syy *= 0.25;

  const double theta = 0.5 * std::atan2( 2.0 * sxy, sxx - syy );
  r.ux               = std::cos( theta );
  r.uy               = std::sin( theta );
  r.vx               = -r.uy;
  r.vy               = r.ux;

  double min_u = std::numeric_limits<double>::infinity();
  double max_u = -std::numeric_limits<double>::infinity();
  double min_v = std::numeric_limits<double>::infinity();
  double max_v = -std::numeric_limits<double>::infinity();

  for( int i = 0; i < 4; ++i )
  {
    const double dx = corners[i].x - r.cx;
    const double dy = corners[i].y - r.cy;

    const double pu = dot2( dx, dy, r.ux, r.uy );
    const double pv = dot2( dx, dy, r.vx, r.vy );

    min_u = std::min( min_u, pu );
    max_u = std::max( max_u, pu );
    min_v = std::min( min_v, pv );
    max_v = std::max( max_v, pv );
  }

  r.hu = std::max( 0.0, 0.5 * ( max_u - min_u ) );
  r.hv = std::max( 0.0, 0.5 * ( max_v - min_v ) );
  return true;
}

inline double
rect_support_along_dir( const OrientedRect& r, double nx, double ny )
{
  return r.hu * std::abs( dot2( r.ux, r.uy, nx, ny ) ) + r.hv * std::abs( dot2( r.vx, r.vy, nx, ny ) );
}

// =============================================================================
// lane band selection
// =============================================================================
struct LaneBandDescriptor
{
  std::shared_ptr<adore::map::Lane> left_lane;
  std::shared_ptr<adore::map::Lane> right_lane;
  bool                              left_use_outer  = true;
  bool                              right_use_outer = true;
};

LaneBandDescriptor
make_lane_band( const std::map<double, std::shared_ptr<adore::map::Lane>>& lane_map, const std::shared_ptr<adore::map::Lane>& route_lane,
                LaneScope lane_scope )
{
  LaneBandDescriptor desc;

  auto fallback_to_my_lane = [&]() {
    desc.left_lane       = route_lane;
    desc.right_lane      = route_lane;
    desc.left_use_outer  = route_lane->left_of_reference;
    desc.right_use_outer = !route_lane->left_of_reference;
  };

  fallback_to_my_lane();

  if( lane_map.empty() || lane_map.size() == 1U || lane_scope == LaneScope::MyLane )
    return desc;

  const bool route_left_side = route_lane->left_of_reference;

  const auto global_min_lane = lane_map.begin()->second;
  const auto global_max_lane = lane_map.rbegin()->second;

  std::shared_ptr<adore::map::Lane> inner_left_lane  = nullptr; // first lane with offset > 0
  std::shared_ptr<adore::map::Lane> inner_right_lane = nullptr; // last lane with offset < 0

  auto first_pos_it = lane_map.upper_bound( 0.0 );
  if( first_pos_it != lane_map.end() )
    inner_left_lane = first_pos_it->second;

  auto first_nonneg_it = lane_map.lower_bound( 0.0 );
  if( first_nonneg_it != lane_map.begin() )
    inner_right_lane = std::prev( first_nonneg_it )->second;

  switch( lane_scope )
  {
    case LaneScope::SameDirection:
    {
      const auto left_boundary_lane  = route_left_side ? global_max_lane : inner_right_lane;
      const auto right_boundary_lane = route_left_side ? inner_left_lane : global_min_lane;

      if( left_boundary_lane && right_boundary_lane )
      {
        desc.left_lane       = left_boundary_lane;
        desc.right_lane      = right_boundary_lane;
        desc.left_use_outer  = route_left_side;
        desc.right_use_outer = !route_left_side;
      }
      break;
    }

    case LaneScope::AllLanes:
    default:
    {
      if( global_max_lane && global_min_lane )
      {
        desc.left_lane       = global_max_lane;
        desc.right_lane      = global_min_lane;
        desc.left_use_outer  = true;
        desc.right_use_outer = true;
      }
      break;
    }
  }

  return desc;
}

// =============================================================================
// scoped timing
// =============================================================================
// using clock     = std::chrono::steady_clock;
// using ms_double = std::chrono::duration<double, std::milli>;

// struct ScopedTimer
// {
//   const char*       label;
//   clock::time_point t0;

// ScopedTimer( const char* l ) :
//   label( l ),
//   t0( clock::now() )
// {}

// ~ScopedTimer()
// {
//   const auto dt = ms_double( clock::now() - t0 ).count();
//   std::cerr << "[DrivableArea] " << label << " took " << dt << " ms\n";
// }
// };

} // anonymous namespace

// =============================================================================
// boundary initialization
// =============================================================================
void
initialize_boundaries( DrivableArea& area, const adore::map::Route& route, double start_s, double end_s, LaneScope lane_scope )
{
  area.left_boundary.clear();
  area.right_boundary.clear();

  auto it  = route.reference_line.lower_bound( start_s );
  auto end = route.reference_line.upper_bound( end_s );

  for( ; it != end; ++it )
  {
    const double                s         = it->first;
    const adore::map::MapPoint& map_point = it->second; // contains lane-relative s (per your model)

    if( !area.left_boundary.empty() && std::abs( s - area.left_boundary.rbegin()->first ) < 0.1 )
      continue; // skip too-close points

    const auto frame = make_route_frame( route, s );

    const auto& mp      = it->second;
    auto        lane_it = route.map->lanes.find( mp.parent_id );
    if( lane_it == route.map->lanes.end() )
      continue;

    auto road_it = route.map->roads.find( lane_it->second->road_id );
    if( road_it == route.map->roads.end() )
      continue;

    const auto band = make_lane_band( road_it->second.lane_offset_to_lane, lane_it->second, lane_scope );

    const auto& left_border  = band.left_use_outer ? band.left_lane->borders.outer : band.left_lane->borders.inner;
    const auto& right_border = band.right_use_outer ? band.right_lane->borders.outer : band.right_lane->borders.inner;

    auto left_mp  = left_border.get_interpolated_point( map_point.s );
    auto right_mp = right_border.get_interpolated_point( map_point.s );

    const double l_left  = lateral_offset( frame, left_mp.x, left_mp.y );
    const double l_right = lateral_offset( frame, right_mp.x, right_mp.y );

    if( l_left >= l_right )
    {
      area.left_boundary[s]  = left_mp;
      area.right_boundary[s] = right_mp;
    }
    else
    {
      area.left_boundary[s]  = right_mp;
      area.right_boundary[s] = left_mp;
    }
  }
}

void
carve_participants( DrivableArea& area, const adore::map::Route& route, double start_s, double end_s,
                    const dynamics::TrafficParticipantSet& traffic_participants, const DrivableAreaConfig& config )
{
  if( area.left_boundary.empty() || area.right_boundary.empty() )
    return;

  for( const auto& [id, participant] : traffic_participants.participants )
  {
    if( participant.state.vx > 0.1 )
      continue;

    auto s_opt = route.get_s( participant.state, k_participant_s_range );
    if( !s_opt )
      continue;

    const double participant_s = *s_opt;
    if( participant_s < start_s || participant_s > end_s )
      continue;

    const auto corners_poly = participant.get_corners( config.longitudinal_inflation, config.lateral_inflation );
    if( corners_poly.points.size() < 4 )
      continue;

    const auto& corners = corners_poly.points;

    // longitudinal span along route (global s)
    double s_min = std::numeric_limits<double>::infinity();
    double s_max = -std::numeric_limits<double>::infinity();
    for( const auto& c : corners )
    {
      auto cs = route.get_s( c, k_participant_s_range );
      if( cs )
      {
        s_min = std::min( s_min, *cs );
        s_max = std::max( s_max, *cs );
      }
    }

    if( !std::isfinite( s_min ) || !std::isfinite( s_max ) )
      continue;

    s_min = std::max( s_min - k_small_long_buffer, start_s );
    s_max = std::min( s_max + k_small_long_buffer, end_s );
    if( s_min > s_max )
      continue;

    OrientedRect rect;
    if( !make_oriented_rect_from_corners( corners, rect ) )
      continue;

    auto it_end = area.left_boundary.upper_bound( s_max );
    for( auto it = area.left_boundary.lower_bound( s_min ); it != it_end; ++it )
    {
      const double s_key = it->first;

      auto right_it = area.right_boundary.find( s_key );
      if( right_it == area.right_boundary.end() )
        continue;

      auto frame = make_route_frame( route, s_key );

      auto& left_mp  = it->second;
      auto& right_mp = right_it->second;

      // compute l of actual points
      double l_left  = lateral_offset( frame, left_mp.x, left_mp.y );
      double l_right = lateral_offset( frame, right_mp.x, right_mp.y );

      // IMPORTANT: keep point ordering consistent (left is max-l, right is min-l)
      if( l_left < l_right )
      {
        std::swap( left_mp, right_mp );
        std::swap( l_left, l_right );
      }

      const double corridor_width = l_left - l_right;
      if( corridor_width < k_min_width )
        continue;

      // participant interval along route normal at this station
      const double center_l = lateral_offset( frame, rect.cx, rect.cy );
      const double radius_l = rect_support_along_dir( rect, frame.nx, frame.ny );

      const double a = center_l - radius_l;
      const double b = center_l + radius_l;

      const double o_min = std::max( a, l_right );
      const double o_max = std::min( b, l_left );
      if( o_min >= o_max )
        continue;

      const bool blocks_all = ( o_min <= l_right + k_station_epsilon ) && ( o_max >= l_left - k_station_epsilon );

      double new_l_right = l_right;
      double new_l_left  = l_left;

      if( blocks_all )
      {
        // participant fully covers corridor at this s -> collapse to 0 width (single line)
        const double mid = 0.5 * ( l_left + l_right );
        new_l_right      = mid;
        new_l_left       = mid;
      }
      else
      {
        const double free_left  = l_left - o_max;  // keep [o_max, l_left]
        const double free_right = o_min - l_right; // keep [l_right, o_min]

        const bool can_keep_left  = free_left > k_station_epsilon;
        const bool can_keep_right = free_right > k_station_epsilon;

        if( can_keep_left && !can_keep_right )
        {
          // only left passage exists -> cut right side (move right boundary inward)
          new_l_right = o_max;
        }
        else if( can_keep_right && !can_keep_left )
        {
          // only right passage exists -> cut left side (move left boundary inward)
          new_l_left = o_min;
        }
        else if( can_keep_left && can_keep_right )
        {
          // participant fully inside corridor -> keep the larger passage
          if( free_left >= free_right )
          {
            new_l_right = o_max; // keep left passage
          }
          else
          {
            new_l_left = o_min; // keep right passage
          }
        }
        else
        {
          // extremely degenerate numeric case: overlap but both free widths ~0
          const double mid = 0.5 * ( l_left + l_right );
          new_l_right      = mid;
          new_l_left       = mid;
        }
      }
    }
  }
}

// =============================================================================
// reference line optimization
// =============================================================================
void
build_reference_line( DrivableArea& area, const adore::map::Route& route, const DrivableAreaConfig& config )
{
  if( area.left_boundary.empty() )
    return;

  struct Station
  {
    double s;
    double l_min;
    double l_max;
    double cx, cy;
    double nx, ny;
    double l_target = 0.0; // target lateral offset in this station’s frame
    bool   blocked  = false;
  };

  std::vector<Station> stations;
  stations.reserve( area.left_boundary.size() );

  for( const auto& [s, left_mp] : area.left_boundary )
  {
    auto it = area.right_boundary.find( s );
    if( it == area.right_boundary.end() )
      continue;

    const auto frame = make_route_frame( route, s );
    if( !frame.ok )
      continue;

    const auto& right_mp = it->second;

    Station st;
    st.s  = s;
    st.cx = frame.cx;
    st.cy = frame.cy;
    st.nx = frame.nx;
    st.ny = frame.ny;

    st.l_max = lateral_offset( frame, left_mp.x, left_mp.y ) - config.safety_margin;
    st.l_min = lateral_offset( frame, right_mp.x, right_mp.y ) + config.safety_margin;
    if( st.l_min > st.l_max )
      std::swap( st.l_min, st.l_max );

    if( math::distance_2d( left_mp, right_mp ) < 2.0 )
    {
      const double mid = 0.5 * ( st.l_min + st.l_max );
      st.l_min = st.l_max = mid;
      st.blocked          = true;
    }

    stations.push_back( st );
  }

  const int n = static_cast<int>( stations.size() );
  if( n < 3 )
    return;

  Eigen::SparseMatrix<double> H( n, n );
  Eigen::VectorXd             g = Eigen::VectorXd::Zero( n );

  constexpr double w_target    = 1.0;
  constexpr double w_rate      = 1.0;
  constexpr double w_curvature = 1.0;

  for( int i = 0; i < n; ++i )
  {
    H.coeffRef( i, i ) += w_target;
    g( i )             += -w_target * stations[i].l_target;
  }


  for( int i = 0; i < n - 1; ++i )
  {
    const double ds = std::max( stations[i + 1].s - stations[i].s, k_station_epsilon );
    const double w  = w_rate / ( ds * ds );

    H.coeffRef( i, i )         += w;
    H.coeffRef( i + 1, i + 1 ) += w;
    H.coeffRef( i, i + 1 )     -= w;
    H.coeffRef( i + 1, i )     -= w;
    H.coeffRef( i, i )         += 1.0;
  }

  for( int i = 1; i < n - 1; ++i )
  {
    const double ds1 = std::max( stations[i].s - stations[i - 1].s, k_station_epsilon );
    const double ds2 = std::max( stations[i + 1].s - stations[i].s, k_station_epsilon );
    const double w   = w_curvature / std::pow( 0.5 * ( ds1 + ds2 ), 4 );

    H.coeffRef( i - 1, i - 1 ) += w;
    H.coeffRef( i, i )         += 4.0 * w;
    H.coeffRef( i + 1, i + 1 ) += w;
    H.coeffRef( i - 1, i )     -= 2.0 * w;
    H.coeffRef( i, i - 1 )     -= 2.0 * w;
    H.coeffRef( i, i + 1 )     -= 2.0 * w;
    H.coeffRef( i + 1, i )     -= 2.0 * w;
    H.coeffRef( i - 1, i + 1 ) += w;
    H.coeffRef( i + 1, i - 1 ) += w;
  }

  Eigen::SparseMatrix<double> A( n, n );
  A.setIdentity();

  Eigen::VectorXd lb( n ), ub( n );
  for( int i = 0; i < n; ++i )
  {
    lb( i ) = stations[i].l_min;
    ub( i ) = stations[i].l_max;
  }

  OsqpEigen::Solver solver;
  solver.settings()->setVerbosity( false );
  solver.settings()->setWarmStart( true );
  solver.data()->setNumberOfVariables( n );
  solver.data()->setNumberOfConstraints( n );
  solver.data()->setHessianMatrix( H );
  solver.data()->setGradient( g );
  solver.data()->setLinearConstraintsMatrix( A );
  solver.data()->setLowerBound( lb );
  solver.data()->setUpperBound( ub );

  if( !solver.initSolver() || solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError )
    return;

  const auto sol = solver.getSolution();

  area.reference_line.clear();
  for( int i = 0; i < n; ++i )
  {
    adore::map::MapPoint mp;
    mp.s = stations[i].s;
    mp.x = stations[i].cx + sol( i ) * stations[i].nx;
    mp.y = stations[i].cy + sol( i ) * stations[i].ny;
    if( stations[i].blocked )
      mp.max_speed = 0.0;

    area.reference_line[mp.s] = mp;
  }
}

// =============================================================================
// public API
// =============================================================================
DrivableArea
create_drivable_area( const adore::map::Route& route, double start_s, double end_s,
                      const dynamics::TrafficParticipantSet& traffic_participants, const DrivableAreaConfig& config )
{

  DrivableArea area;
  if( !route.map )
    return area;

  if( end_s < start_s )
    std::swap( start_s, end_s );

  {
    // ScopedTimer t( "initialize_boundaries" );
    initialize_boundaries( area, route, start_s, end_s, config.lane_scope );
  }

  if( area.left_boundary.empty() )
    return area;

  {
    // ScopedTimer t( "carve_participants" );
    carve_participants( area, route, start_s, end_s, traffic_participants, config );
  }

  {
    // ScopedTimer t( "build_reference_line" );
    build_reference_line( area, route, config );
  }

  return area;
}

} // namespace planner
} // namespace adore
