% Projects a symmetric matrix on to the set of positive definite matrices with
% rank less than or equal to k.

function X = project_psd_rank(X, k)
  [V, D] = eig(X);
  d = diag(D);
  d = max(0, d);
  [d, i] = sort(d, 'descend');
  V = V(:, i);
  if length(d) > 3
    d(4:end) = 0;
  end
  X = V * diag(d) * V';
end
