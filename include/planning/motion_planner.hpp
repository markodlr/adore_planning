#pragma once

#include <map>
#include <memory>
#include <optional>
#include <string>

#include "adore_map/route.hpp"
#include "dynamics/comfort_settings.hpp"
#include "dynamics/traffic_participant.hpp"
#include "dynamics/vehicle_state.hpp"
#include "planning/common/planning_config.hpp"
#include "planning/drivable_area.hpp"
#include "planning/nominal_participant_prediction.hpp"
#include "planning/trajectory_optimizer.hpp"

namespace adore
{
namespace planner
{

struct MotionPlannerConfig
{
  DrivableAreaConfig da_config;
  SpeedProfileConfig sp_config;
  
  // General planning settings can be added here
  double drivable_area_length = 100.0;
  double drivable_area_before = 10.0;
};

struct PlannerResult
{
  std::optional<dynamics::Trajectory> trajectory;
  std::optional<DrivableArea>         drivable_area; // For debugging/visualization
  bool                                success = false;
  std::string                         message;
};

class MotionPlanner
{
public:
  MotionPlanner();

  /*
   * Initialize the planner with configuration and vehicle parameters.
   */
  void init( const MotionPlannerConfig& config, const dynamics::PhysicalVehicleParameters& vehicle_params,
             const std::shared_ptr<dynamics::ComfortSettings>& comfort_settings );

  /*
   * Update internal parameters from a key-value map (e.g. from ROS parameters).
   * This is useful for tuning parameters at runtime or load time.
   */
  void update_parameters( const std::map<std::string, double>& params );

  /*
   * Main planning function.
   * Generates a trajectory to follow the given route while avoiding obstacles.
   */
  PlannerResult plan( const adore::map::Route& route, const dynamics::VehicleStateDynamic& ego_state,
                      const dynamics::TrafficParticipantSet& participants );
                      
  MotionPlannerConfig                        config;

private:
  dynamics::PhysicalVehicleParameters        vehicle_params;
  std::shared_ptr<dynamics::ComfortSettings> comfort_settings;

  // Sub-components
  TrajectoryOptimizer          optimizer;
  NominalParticipantPrediction participant_predictor;
};

} // namespace planner
} // namespace adore
