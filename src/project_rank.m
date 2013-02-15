% Projects a matrix on to the set with rank <= k.

function X = project_rank(X, k)
  [U, S, V] = svd(X, 'econ');

  U = U(:, 1:k);
  V = V(:, 1:k);
  s = diag(S);
  s = s(1:k);
  S = diag(s);

  X = U * S * V';
end
