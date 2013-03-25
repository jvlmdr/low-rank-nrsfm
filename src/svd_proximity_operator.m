% Solves
% arg min_{M, B} 1/2 ||W - M B||_F^2 + 1/2 rho ||M - V||_F^2
%
% W -- m x n
% V -- m x r (r is rank of factorization)

function [M, B] = svd_proximity_operator(W, V, rho)
  m = size(W, 1);
  n = size(W, 2);
  r = size(V, 2);

  % Remove rho.
  V = sqrt(rho) * V;

  X = [W, V];
  Y = X * X';
  [U, D] = eig(Y);

  % Sort eigenvalues.
  d = diag(D);
  [~, ind] = sort(abs(d), 'descend');
  d = d(ind);
  U = U(:, ind);

  % Take largest.
  Q = U(:, 1:r);

  % Extract solution.
  B_hat = Q' * W;
  A = Q' * V;
  B = inv(A) * B_hat;
  M = Q * A;

  % Reinstate rho.
  M = 1 / sqrt(rho) * M;
  B = sqrt(rho) * B;
end
