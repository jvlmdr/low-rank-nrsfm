function X = k_reshape(X, k)
  [km, n] = size(X);
  m = km / k;
  % [km, n] -> [k, m, n]
  X = reshape(X, [k, m, n]);
  % [k, m, n] -> [k, n, m]
  X = permute(X, [1, 3, 2]);
  % [k, n, m] -> [kn, m]
  X = reshape(X, [k * n, m]);
end
