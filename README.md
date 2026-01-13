# Planning Library for Autonomous Vehicles

## Overview
The **Planning Library** provides core algorithms and tools for trajectory planning, speed profiling, and optimization in autonomous vehicles. It is designed to handle complex scenarios including obstacle avoidance, lane following, and multi-agent coordination.

---

## Features
- **Drivable Area Generation**:
  - Computes safe corridors (boundaries) and reference lines based on the map route and static/dynamic obstacles.
- **Speed Profile Planning**:
  - Generates longitudinal speed profiles ($s(t)$) using ST-graph dynamic programming followed by quadratic programming smoothing.
- **Trajectory Optimization**:
  - Optimizes vehicle trajectories using Model Predictive Control (MPC) to ensure dynamic feasibility and comfort.
- **Multi-Agent Planning**:
  - Coordinates multiple traffic participants using game-theoretic or cooperative planning approaches.
- **Prediction**:
  - Provides nominal predictions for other traffic participants.

---

## Included Modules

### Drivable Area
**File:** `include/planning/drivable_area.hpp`
- Constructs a `DrivableArea` object containing left/right boundaries and a reference line.
- Accounts for the static route geometry and dynamic traffic participants to define the safe maneuverable space.
- Supports candidate generation for multiple topological choices (e.g., nudging left vs. right).

### Speed Profile
**File:** `include/planning/speed_profile.hpp`
- Plans a longitudinal velocity profile along the `DrivableArea`'s reference line.
- **Methodology**:
  1.  **ST-Graph Construction**: Builds a Space-Time graph considering obstacles.
  2.  **Dynamic Programming**: Finds a coarse optimal path through the discretized ST-graph.
  3.  **Convex Tunnel**: Constructs a collision-free convex corridor around the DP path.
  4.  **QP Smoothing**: Optimizes the final speed profile ($s, v, a, j$) within the tunnel for smoothness and comfort.

### Trajectory Optimizer
**File:** `include/planning/trajectory_optimizer.hpp`
- A single-agent trajectory optimizer based on MPC.
- Minimizes cost functions related to lane deviation, speed tracking, heading error, and control efforts (steering, acceleration).
- Uses a kinematic bicycle model for vehicle dynamics.

### Multi-Agent Planner
**File:** `include/planning/multi_agent_planner.hpp`
- Handles planning for multiple agents simultaneously.
- Can solve for cooperative behaviors or interaction-aware planning.
- Utilizes the underlying `multi_agent_solver` framework.

### Prediction & Control
- **Nominal Participant Prediction** (`nominal_participant_prediction.hpp`): specific logic for predicting future states of non-ego participants.
- **Multi-Agent PID** (`multi_agent_PID.hpp`): A PID-based controller approach for multi-agent scenarios
