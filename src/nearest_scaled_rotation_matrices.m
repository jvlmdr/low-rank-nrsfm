% [R, c] = nearest_scaled_rotation_matrices(M)
%
% Parameters:
% M -- 2F x 3
%
% Returns:
% R -- 2 x 3 x F
% c -- F x 1

function [R, c] = nearest_scaled_rotation_matrices(M)

  F = size(M, 1) / 2;

  % [2F, 3] -> [2, F, 3] -> [2, 3, F]
  M = permute(reshape(M, [2, F, 3]), [1, 3, 2]);

  R = zeros(2, 3, F);
  c = zeros(F, 1);

  for t = 1:F
    [R(:, :, t), c(t)] = nearest_scaled_rotation_matrix(M(:, :, t));
  end
end
