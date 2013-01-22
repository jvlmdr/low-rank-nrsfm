function [A, b] = projection_equations_to_single_equation(projections)
  num_frames = projections.num_frames;

  trajectories = projection_equations_to_trajectory_equations(projections);
  [A, b] = trajectory_equations_to_single_equation(trajectories, num_frames);
end
