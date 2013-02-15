function X = k_reshape(X, k)
  [km, n] = size(X);
  m = km / k;
  X = reshape(X, [k, m, n]);
  X = permute(X, [1, 3, 2]);
  X = reshape(X, [k * n, m]);
end
