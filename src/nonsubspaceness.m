% Parameters:
% X -- 3 x n

function y = nonsubspaceness(X)
  [U, S, V] = svd(X, 'econ');
  A = U(:, 1:2);
  y = norm(X - A * (A' * X), 'fro') / norm(X, 'fro');
end
