% Parameters:
% points -- num_frames x num_points x num_dims

function points = subtract_mean(points)
  num_frames = size(points, 1);
  num_points = size(points, 2);
  num_dims = size(points, 3);

  num_total = num_frames * num_points;

  % Subtract mean.
  points = reshape(points, [num_total, num_dims]);
  mu = mean(points)';
  % Subtract mean.
  points = points - ones(num_total, 1) * mu';
  % Restore dimension.
  points = reshape(points, [num_frames, num_points, num_dims]);
end
