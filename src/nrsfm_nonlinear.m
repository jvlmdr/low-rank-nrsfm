% Parameters:
% W -- 2 x P x F
% R_init -- 2 x 3 x F
% S_init -- 3 x P x F
% K -- Rank of structure
%
% Returns:
% R -- 2 x 3 x F
% S -- 3 x P x F

function [R, S, B, C] = nrsfm_nonlinear(W, R_init, S_init, K, max_iter, tol)
  assert(ndims(W) == 3);
  assert(size(W, 1) == 2);
  P = size(W, 2);
  F = size(W, 3);

  assert(ndims(R_init) == 3);
  assert(all(size(R_init) == [2, 3, F]));

  assert(ndims(S_init) == 3);
  assert(all(size(S_init) == [3, P, F]));

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

  % Compute initial residual for debug purposes. Should reassure us that
  % quaternion conversion worked and matrices are being addressed correctly.
  W_mat = reshape(permute(W, [1, 3, 2]), [2 * F, P]);
  R_mat = block_diagonal_cameras(R_init);
  B_mat = reshape(permute(B_init, [1, 3, 2]), [3 * P, K]);
  S_mat = B_mat * C_init;
  S_mat = reshape(permute(reshape(S_mat, [3, P, F]), [1, 3, 2]), [3 * F, P]);
  fprintf('Initial residual: %g\n', ...
      1/2 * norm(W_mat - R_mat * S_mat, 'fro') ^ 2);

  % Solve.
  [Q, B, C] = nrsfm_nonlinear_mex(W, Q_init, B_init, C_init, max_iter, tol);

  % Compose structure.
  B_mat = reshape(permute(B, [1, 3, 2]), [3 * P, K]);
  S = B_mat * C;
  S = reshape(S, [3, P, F]);

  % Convert back to rotation matrices.
  R = zeros(2, 3, F);
  for t = 1:F
    Rt = quat2rot(Q(:, t));
    R(:, :, t) = Rt(1:2, :);
  end

  % Compute final residual for debug purposes.
  R_mat = block_diagonal_cameras(R);
  S_mat = reshape(permute(S, [1, 3, 2]), [3 * F, P]);
  fprintf('Final residual: %g\n', 1/2 * norm(W_mat - R_mat * S_mat, 'fro') ^ 2);
end
