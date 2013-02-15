function [R, c] = nearest_scaled_rotation_matrices(P)

  F = size(P, 1) / 2;

  P = permute(reshape(P, [2, F, 3]), [1, 3, 2]);

  R = zeros(2, 3, F);
  c = zeros(F, 1);

  for t = 1:F
    [R(:, :, t), c(t)] = nearest_scaled_rotation_matrix(P(:, :, t));
  end
end
