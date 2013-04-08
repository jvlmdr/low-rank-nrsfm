% Parameters:
% projections -- 2 x P x F
% structure -- 3 x P x F
%
% Returns:
% rotations -- 2 x 3 x F

function rotations = find_cameras(projections, structure)
  F = size(projections, 3);
  rotations = zeros(2, 3, F);

  for t = 1:F
    W_t = projections(:, :, t);
    S_t = structure(:, :, t);

    Q = W_t / S_t;
    Q = nearest_scaled_rotation_matrix(Q);
    rotations(:, :, t) = Q(1:2, :);
  end
end
