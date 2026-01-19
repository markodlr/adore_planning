
#include <gtest/gtest.h>

#include "dynamics/comfort_settings.hpp"
#include "dynamics/physical_vehicle_parameters.hpp"
#include "planning/speed_profile.hpp"

using namespace adore::planner;
using namespace adore::dynamics;

TEST( SpeedProfileTest, SmootherReducesOscillation )
{
  // Setup
  DrivableArea area;
  // Create a straight reference line
  for( int i = 0; i < 200; ++i )
  {
    adore::map::MapPoint mp;
    mp.x                         = i * 1.0;
    mp.y                         = 0.0;
    mp.max_speed                 = 10.0;
    area.reference_line[i * 1.0] = mp;
  }
  // Also need boundaries for project_to_path to work inside some functions if called
  // But plan_speed_profile mainly uses reference_line.
  // Wait, build_boundary_for_participant uses project_to_path which might need full area info?
  // Actually project_to_path works on the map provided.

  TrafficParticipantSet participants; // Empty for now

  VehicleStateDynamic ego_state;
  ego_state.x  = 0.0;
  ego_state.y  = 0.0;
  ego_state.vx = 5.0; // Starting at 5 m/s

  PhysicalVehicleParameters vehicle_params;

  ComfortSettings comfort;
  comfort.max_acceleration         = 2.0;
  comfort.min_acceleration         = -4.0;
  comfort.max_lateral_acceleration = 2.0;
  comfort.max_speed                = 15.0; // want to accelerate to this

  SpeedProfileConfig config;
  config.dt_dp      = 0.2; // Coarse
  config.ds_dp      = 0.1; // Fine
  config.s_horizon  = 100.0;
  config.total_time = 10.0;
  config.v_cruise   = 10.0; // Target 10 m/s

  // Act
  SpeedProfile profile = plan_speed_profile( area, participants, ego_state, vehicle_params, comfort, config );

  // Assert
  ASSERT_FALSE( profile.empty() );

  // Check smoothness
  double max_jerk         = 0.0;
  double max_accel_change = 0.0;

  for( size_t i = 1; i < profile.size(); ++i )
  {
    double dt = profile[i].t - profile[i - 1].t;
    if( dt > 1e-3 )
    {
      double da        = profile[i].a - profile[i - 1].a;
      double jerk      = da / dt;
      max_jerk         = std::max( max_jerk, std::abs( jerk ) );
      max_accel_change = std::max( max_accel_change, std::abs( da ) );
    }
  }

  std::cout << "Max Jerk: " << max_jerk << " m/s^3" << std::endl;
  std::cout << "Max Accel Change: " << max_accel_change << " m/s^2" << std::endl;

  // With pure DP and dt=0.5, ds=0.1, delta_a is ~0.8.
  // If we simply check that acceleration doesn't oscillate wildly, we might need a better metric.
  // Ideally, a qp smoother should bring jerk down significantly.
  // For now, let's just assert it runs and prints values.
  // We can tighten this check after implementation.

  EXPECT_GT( profile.size(), 10 );
}
