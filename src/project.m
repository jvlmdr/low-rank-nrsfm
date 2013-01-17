% Projects 3D points to 2D points.
%
% Parameters:
% points -- num_frames x num_points x 3

function projections = project(cameras, points)
  num_frames = size(points, 1);
  num_points = size(points, 2);

  for i = 1:num_frames
    % Get points in this frame.
    % size(X) => [3, num_points]
    X = shiftdim(points(i, :, :), 1)';

    % Make homogeneous.
    X = [X; ones(1, num_points)];
    % Project.
    W = cameras(:, :, i) * X;
    % Make inhomogeneous.
    W = W(1:2, :) * spdiags(1 ./ W(3, :)', 0, num_points, num_points);

    projections(i, :, :) = W';
  end
end
