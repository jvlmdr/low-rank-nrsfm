function S = reconstruct_nullspace(W, K)
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

  % Solve for rotations using a single triple.
  Rs = find_rotations(M_hat, 1e6);

  R = block_diagonal_cameras(Rs);
  [G, C] = find_corrective_matrix_nullspace(M_hat, Rs);
  S = kron(C, eye(3)) * inv(G) * B_hat;

  % [3F, P] -> [3, F, P] -> [F, P, 3]
  S = reshape(S, [3, F, P]);
  S = permute(S, [2, 3, 1]);
end
