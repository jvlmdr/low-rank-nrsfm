function plot_trajectories(fig, trajectories, num_frames, fps, lim, opts)

  num_points = length(trajectories);
  d = size(trajectories(1).points, 2);

  points = nan(num_frames, num_points, d);

  for i = 1:num_points
    points(trajectories(i).frames, i, :) = trajectories(i).points;
  end

  plot_movie(fig, points, fps, lim, opts);
end
