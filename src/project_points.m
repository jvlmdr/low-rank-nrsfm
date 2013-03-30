% Takes a series of 3D points and cameras and does the projection.
%
% Parameters:
% cameras -- F-length vector of camera structs.
% points -- 3 x N x F matrix of 3D points.
%   Or 3 x N matrix of 3D points for rigid case.
%
% Returns:
% 2 x N x F

function projections = project_points(cameras, points)
  F = length(cameras);
  N = size(points, 2);

  % Initialise projection results.
  projections = zeros(2, N, F);

  for i = 1:F
    P = cameras(i).P;
    if length(size(points)) == 3
      projections(:, :, i) = project(P, points(:, :, i));
    elseif length(size(points)) == 2
      projections(:, :, i) = project(P, points);
    end
  end
end
