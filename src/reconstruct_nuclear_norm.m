function S = reconstruct_nuclear_norm(W, K)
  F = size(W, 1);
  P = size(W, 2);

  % [F, P, 2] -> [2, F, P] -> [2F, P]
  W = permute(W, [3, 1, 2]);
  W = reshape(W, [2 * F, P]);

  % Subtract centroid per frame.
  mu = mean(W, 2);
  W = W - mu * ones(1, P);

  % Compute SVD of W.
  [U, D, V] = svd(W, 'econ');
  U = U(:, 1:3 * K);
  V = V(:, 1:3 * K);
  d = diag(D);
  d = d(1:3 * K);

  % Get initial factorization.
  M_hat = 1 / sqrt(d(1)) * U * diag(sqrt(d));
  B_hat = sqrt(d(1)) * diag(sqrt(d)) * V';

  Rs = find_rotations(M_hat, 1e6);
  S = find_structure_affine_cameras(W, Rs, true, ...
      struct(...
        'rho', 1, ...
        'mu', 10, ...
        'tau_incr', 2, ...
        'tau_decr', 2, ...
        'max_iter', 80, ...
        'epsilon_abs', 1e-3, ...
        'epsilon_rel', 1e-3, ...
        'min_rho_iter', 4));

  % [3F, P] -> [3, F, P] -> [F, P, 3]
  S = reshape(S, [3, F, P]);
  S = permute(S, [2, 3, 1]);
end
