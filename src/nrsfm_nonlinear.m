% Parameters:
% W -- 2 x P x F
% R_init -- 2 x 3 x F
% S_init -- 3 x P x F
% K -- Rank of structure
%
% Returns:
% R -- 2 x 3 x F
% S -- 3 x P x F

function [R, S] = nrsfm_nonlinear(W, R_init, S_init, K, max_iter, tol)
  P = size(W, 2);
  F = size(W, 3);

  % Convert rotations to quaternions.
  Q_init = zeros(4, F);
  for t = 1:F
    Rt = R_init(:, :, t);
    Rt = [Rt; cross(Rt(1, :), Rt(2, :))];
    Q_init(:, t) = rot2quat(Rt);
  end

  % Convert structure to basis and coefficients.
  S = reshape(S_init, [3 * P, F]);
  [U, D, V] = svd(S);
  % 3 x K x P basis vectors
  B_init = U(:, 1:K);
  B_init = permute(reshape(B_init, [3, P, K]), [1, 3, 2]);
  % K x F coefficients
  C_init = D(1:K, 1:K) * V(:, 1:K)';

  % Solve.
  [Q, B, C] = nrsfm_nonlinear_mex(W, Q_init, B_init, C_init, max_iter, tol);

  % Compose structure.
  B = reshape(permute(B, [1, 3, 2]), [3 * P, K]);
  S = B * C;
  S = reshape(S, [3, P, F]);

  % Convert back to rotation matrices.
  R = zeros(2, 3, F);
  for t = 1:F
    Rt = quat2rot(Q(:, t));
    R(:, :, t) = Rt(1:2, :);
  end
end
