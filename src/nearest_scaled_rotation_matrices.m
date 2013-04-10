% [R, c] = nearest_scaled_rotation_matrices(M)
%
% Parameters:
% Q -- 2 x 3 x F
%
% Returns:
% R -- 2 x 3 x F
% c -- F x 1

function [R, c] = nearest_scaled_rotation_matrices(Q)
  F = size(Q, 3);
  R = zeros(2, 3, F);
  c = zeros(F, 1);

  for t = 1:F
    [R(:, :, t), c(t)] = nearest_scaled_rotation_matrix(Q(:, :, t));
  end
end
