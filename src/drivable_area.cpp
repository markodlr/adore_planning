/********************************************************************************
 * Copyright (c) 2025 Contributors
 *
 * SPDX-License-Identifier: EPL-2.0
 ********************************************************************************/

#include "planning/drivable_area.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <OsqpEigen/OsqpEigen.h>

#include "adore_map/map.hpp" // For Lane (needed for lane utils)
#include "adore_math/geometry/projection.hpp"
#include "adore_map/route.hpp"
#include "dynamics/traffic_participant.hpp"
#include "planning/common/planning_config.hpp"


namespace adore
{
namespace planner
{

namespace
{

constexpr double k_station_epsilon     = 1e-3; // [m]
constexpr double k_min_width           = 0.2;  // [m] 
constexpr double k_participant_s_range = 10.0;

// =============================================================================
// Internal: RouteFrame (formerly route_frame.hpp)
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
  // Small lookahead for heading
  const double k_heading_probe_ds = 0.1;
  const auto p2 = route.get_map_point_at_s( s + k_heading_probe_ds );

  const double dx  = p2.x - p.x;
  const double dy  = p2.y - p.y;
  
  const double angle = std::atan2( dy, dx );
  return { p.x, p.y, -std::sin( angle ), std::cos( angle ), true };
}

inline double
lateral_offset( const RouteFrame& f, double x, double y )
{
  return ( x - f.cx ) * f.nx + ( y - f.cy ) * f.ny;
}

// =============================================================================
// Internal: Lane Utils
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
// Internal: Obstacle Logic
// =============================================================================

std::map<double, std::pair<double, double>>
compute_obstacle_envelope( const adore::map::Route& route, double start_s, double end_s,
                           const dynamics::TrafficParticipant& participant,
                           const dynamics::PhysicalVehicleParameters& vehicle_params,
                           const DrivableAreaConfig&           config )
{
  std::map<double, std::pair<double, double>> envelope;

  // longitudinal span along route (global s)
  double s_min = std::numeric_limits<double>::infinity();
  double s_max = -std::numeric_limits<double>::infinity();
  for( const auto& c : participant.get_corners().points )
  {
    auto cs = route.get_s( c, k_participant_s_range );
    if( cs )
    {
      s_min = std::min( s_min, *cs );
      s_max = std::max( s_max, *cs );
    }
  }

  if( !std::isfinite( s_min ) || !std::isfinite( s_max ) )
    return envelope;

  // Check if the participant overlaps with the requested station range [start_s, end_s]
  if ( s_max < start_s || s_min > end_s )
    return envelope;

  const auto corners_poly = participant.get_corners( config.longitudinal_inflation, config.lateral_inflation );
  if( corners_poly.points.size() < 4 )
    return envelope;

  const auto& corners = corners_poly.points;

  // longitudinal span along route (global s)
  struct SlCorner { double s; double l; };
  std::vector<SlCorner> sl_corners;
  sl_corners.reserve(corners.size());

  s_min = std::numeric_limits<double>::infinity();
  s_max = -std::numeric_limits<double>::infinity();

  for( const auto& c : corners )
  {
    auto cs = route.get_s( c, k_participant_s_range );
    if( cs )
    {
      const auto frame = make_route_frame( route, *cs );
      const double l = lateral_offset( frame, c.x, c.y );
      sl_corners.push_back({ *cs, l });
      s_min = std::min( s_min, *cs );
      s_max = std::max( s_max, *cs );
    }
  }

  if( sl_corners.size() < 3 || !std::isfinite( s_min ) || !std::isfinite( s_max ) )
    return envelope;

  s_min = std::max( s_min, start_s );
  s_max = std::min( s_max, end_s );
  if( s_min > s_max )
    return envelope;

  // Sample at route resolution
  auto it_start = route.reference_line.lower_bound( s_min );
  auto it_end   = route.reference_line.upper_bound( s_max );

  for( auto it = it_start; it != it_end; ++it )
  {
    const double s_query = it->first;
    double o_min = std::numeric_limits<double>::infinity();
    double o_max = -std::numeric_limits<double>::infinity();
    bool found = false;

    // Check intersections of line s = s_query with edges (sl_i -> sl_j)
    for (size_t i = 0; i < sl_corners.size(); ++i) {
        const auto& a = sl_corners[i];
        const auto& b = sl_corners[(i + 1) % sl_corners.size()];
        
        if ((a.s <= s_query && b.s > s_query) || (b.s <= s_query && a.s > s_query)) {
            const double u = (s_query - a.s) / (b.s - a.s);
            const double l_interp = a.l + u * (b.l - a.l);
            o_min = std::min(o_min, l_interp);
            o_max = std::max(o_max, l_interp);
            found = true;
        }
    }

    if (found) {
        envelope[s_query] = { o_min, o_max };
    }
  }
  return envelope;
}

} // namespace

// =============================================================================
// boundary initialization
// =============================================================================
void
initialize_boundaries( DrivableArea& area, const adore::map::Route& route, double start_s, double end_s, const DrivableAreaConfig& config )
{
  area.left_boundary.clear();
  area.right_boundary.clear();

  auto it  = route.reference_line.lower_bound( start_s );
  auto end = route.reference_line.upper_bound( end_s );

  for( ; it != end; ++it )
  {
    const double                s         = it->first;
    const adore::map::MapPoint& map_point = it->second;

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

    // Use shared logic
    const auto band = make_lane_band( road_it->second.lane_offset_to_lane, lane_it->second, config.lane_scope );

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

// =============================================================================
// Public API Implementations
// =============================================================================

void
get_obstacle_envelope( const adore::map::Route&                route,
                       double                                  start_s,
                       double                                  end_s,
                       const dynamics::TrafficParticipant&     participant,
                       const dynamics::PhysicalVehicleParameters& vehicle_params,
                       const DrivableAreaConfig&               config,
                       std::map<double, std::pair<double, double>>& envelope )
{
    envelope = compute_obstacle_envelope(route, start_s, end_s, participant, vehicle_params, config);
}

std::vector<DrivableArea>
generate_candidates( const DrivableArea&                    initial_area,
                     const adore::map::Route&               route,
                     double                                 start_s,
                     double                                 end_s,
                     const dynamics::TrafficParticipantSet& traffic_participants,
                     const dynamics::PhysicalVehicleParameters& vehicle_params,
                     const DrivableAreaConfig&              config )
{
  std::vector<DrivableArea> candidates;
  candidates.push_back( initial_area );

  struct SortedParticipant {
    int id;
    double s;
  };
  std::vector<SortedParticipant> sorted_parts;
  for( const auto& [id, p] : traffic_participants.participants ) {
    auto s_opt = route.get_s( p.state, k_participant_s_range );
    if (s_opt) {
        SortedParticipant sp;
        sp.id = id;
        sp.s = *s_opt;
        sorted_parts.push_back(sp);
    }
  }
  std::sort(sorted_parts.begin(), sorted_parts.end(), [](const SortedParticipant& a, const SortedParticipant& b){ return a.s < b.s; });

  for( const auto& sp : sorted_parts )
  {
    const auto& participant = traffic_participants.participants.at(sp.id);

    auto envelope = compute_obstacle_envelope( route, start_s, end_s, participant, vehicle_params, config );

    if( envelope.empty() )
      continue;

    std::vector<DrivableArea> next_generation;
    next_generation.reserve( candidates.size() * 2 );

    for( const auto& cand : candidates )
    {
      bool interacts      = false;
      bool possible_left  = true;
      bool possible_right = true;
      bool blocks_any     = false;
      double s_blocked    = 0.0;

      for( const auto& [s, obs_range] : envelope )
      {
        auto lit = cand.left_boundary.find( s );
        auto rit = cand.right_boundary.find( s );
        if( lit == cand.left_boundary.end() || rit == cand.right_boundary.end() )
          continue;

        interacts = true;

        const auto frame = make_route_frame( route, s );
        const double l_left_bound  = lateral_offset( frame, lit->second.x, lit->second.y );
        const double l_right_bound = lateral_offset( frame, rit->second.x, rit->second.y );
        
        const double l_max = std::max(l_left_bound, l_right_bound);
        const double l_min = std::min(l_left_bound, l_right_bound);

        const double o_min = std::max( obs_range.first, l_min );
        const double o_max = std::min( obs_range.second, l_max );

        if (o_min >= o_max) continue;

        blocks_any = true;

        const double free_left  = l_max - o_max;
        const double free_right = o_min - l_min;

        // A gap is passable if it's wider than twice the safety margin
        // (which is defined as half_vehicle_width + buffer)
        const double safety_margin = ( vehicle_params.body_width * 0.5 ) + 0.1;
        const double req = 2.0 * safety_margin;

        if ( free_left < req ) {
            possible_left = false;
            s_blocked = s;
        }
        if ( free_right < req ) {
            possible_right = false;
            s_blocked = s;
        }

      }

      if ( !interacts || !blocks_any ) {
        next_generation.push_back( cand );
        continue;
      }

      if ( !possible_left && !possible_right ) {
          continue;
      }

      if ( possible_left ) {
          DrivableArea new_cand = cand;
          for( const auto& [s, obs_range] : envelope ) {
             auto rit = new_cand.right_boundary.find(s);
             if (rit == new_cand.right_boundary.end()) continue;

             auto frame = make_route_frame( route, s );
             double new_l = obs_range.second;
             double old_l = lateral_offset(frame, rit->second.x, rit->second.y);
             
             if ( new_l > old_l ) {
                adore::map::MapPoint mp = rit->second;
                mp.x = frame.cx + new_l * frame.nx;
                mp.y = frame.cy + new_l * frame.ny;
                rit->second = mp;
             }
          }
          next_generation.push_back( new_cand );
      }

      if ( possible_right ) {
          DrivableArea new_cand = cand;
          for( const auto& [s, obs_range] : envelope ) {
             auto lit = new_cand.left_boundary.find(s);
             if (lit == new_cand.left_boundary.end()) continue;

             auto frame = make_route_frame( route, s );
             double new_l = obs_range.first;
             double old_l = lateral_offset(frame, lit->second.x, lit->second.y);
             
             if ( new_l < old_l ) {
                adore::map::MapPoint mp = lit->second; 
                mp.x = frame.cx + new_l * frame.nx;
                mp.y = frame.cy + new_l * frame.ny;
                lit->second = mp;
             }
          }
          next_generation.push_back( new_cand );
      }
    }
    candidates = next_generation;
  }
  
  if (candidates.empty()) {
      candidates.push_back(initial_area);
  }

  return candidates;
}

double score_candidate( const DrivableArea& area, const adore::map::Route& route )
{
    // Simple cost function:
    // 1. Width Cost (maximize width -> minimize -width)
    // 2. Offset Cost (minimize deviation from road center/ref line)
    
    double total_cost = 0.0;
    double count = 0.0;

    for (const auto& [s, l_mp] : area.left_boundary) {
        auto rit = area.right_boundary.find(s);
        if (rit == area.right_boundary.end()) continue;
        const auto& r_mp = rit->second;

        auto frame = make_route_frame( route, s );
        double l_left = lateral_offset(frame, l_mp.x, l_mp.y);
        double l_right = lateral_offset(frame, r_mp.x, r_mp.y);
        
        double width = l_left - l_right;
        double center = (l_left + l_right) * 0.5;

        // Cost terms
        total_cost += -width;           // Prefer wider
        total_cost += std::abs(center); // Prefer centered on route Ref Line
        
        count += 1.0;
    }
    
    return (count > 0) ? (total_cost / count) : 1e9;
}

// =============================================================================
// Reference line smoothing (OSQP)
// =============================================================================
void
build_reference_line( DrivableArea& area, const adore::map::Route& route,
                      const dynamics::PhysicalVehicleParameters& vehicle_params,
                      const DrivableAreaConfig& config )
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

    const double safety_margin = ( vehicle_params.body_width * 0.5 ) + 0.1;
    st.l_max = lateral_offset( frame, left_mp.x, left_mp.y ) - safety_margin;
    st.l_min = lateral_offset( frame, right_mp.x, right_mp.y ) + safety_margin;
    if( st.l_min > st.l_max )
      std::swap( st.l_min, st.l_max );

    const double req = 2.0 * safety_margin;
    const double dist = math::distance_2d( left_mp, right_mp );
    if( dist < req )
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

} // namespace planner
} // namespace adore

namespace adore
{
namespace planner
{

DrivableArea
create_drivable_area( const adore::map::Route& route, double start_s, double end_s,
                      const dynamics::TrafficParticipantSet& traffic_participants,
                      const dynamics::PhysicalVehicleParameters& vehicle_params,
                      const DrivableAreaConfig& config )
{
  DrivableArea area;
  initialize_boundaries( area, route, start_s, end_s, config );

  std::vector<DrivableArea> candidates = generate_candidates( area, route, start_s, end_s, traffic_participants, vehicle_params, config );

  double       best_score = std::numeric_limits<double>::infinity();
  DrivableArea best_area  = area;

  for( const auto& cand : candidates )
  {
    double score = score_candidate( cand, route );
    if( score < best_score )
    {
      best_score = score;
      best_area  = cand;
    }
  }

  area = best_area;
  build_reference_line( area, route, vehicle_params, config );

  return area;
}

} // namespace planner
} // namespace adore

