% Takes a series of 3D points and cameras and does the projection.
%
% Parameters:
% cameras -- F-length vector of camera structs.
% points -- FxNx3 matrix of 3D points.
%   Or Nx3 matrix of 3D points for rigid case.
%
% Returns:
% F x N x 2

function projections = project_points(cameras, points)
  if length(size(points)) == 3
    % [F, N, 3] -> [N, 3, F] for easy access.
    points = shiftdim(points, 1);
  end

  F = length(cameras);
  N = size(points, 1);

  % Initialise projection results.
  projections = zeros(F, N, 2);

  for i = 1:F
    P = cameras(i).P;
    if length(size(points)) == 3
      projections(i, :, :) = project(P, points(:, :, i)')';
    elseif length(size(points)) == 2
      projections(i, :, :) = project(P, points')';
    end
  end
end
