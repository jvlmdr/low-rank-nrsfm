function [M_hat, B_hat, W] = factorize(projections, K)
  P = size(projections, 2);
  F = size(projections, 3);

  % [2, P, F] -> [2, F, P] -> [2F, P]
  W = projections;
  W = permute(W, [1, 3, 2]);
  W = reshape(W, [2 * F, P]);

  % Subtract centroid per frame.
  mu = mean(W, 2);
  W = W - mu * ones(1, P);

  % Compute SVD of W.
  [U, D, V] = svd(W, 'econ');
  d = diag(D);
  U = U(:, 1:(3 * K));
  V = V(:, 1:(3 * K));
  d = d(1:(3 * K));

  % Get initial factorization.
  M_hat = 1 / sqrt(d(1)) * U * diag(sqrt(d));
  B_hat = sqrt(d(1)) * diag(sqrt(d)) * V';
end
