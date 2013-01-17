function plot_cameras(fig, cameras, axis_length)

  ax = gca(fig);

  F = size(cameras, 3);

  % Get camera centers.
  centers = zeros(F, 3);
  for t = 1:F
    P = cameras(:, :, t);
    R = P(:, 1:3);
    d = P(:, 4);
    centers(t, :) = (-R' * d)';
  end

  lim = axis_limits(centers);
  init3(fig, lim);

  for t = 1:F
    plot_camera(ax, cameras(:, :, t), axis_length);
  end
end
