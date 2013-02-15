function X = k_unreshape(X, k)
  [kn, m] = size(X);
  n = kn / k;
  X = reshape(X, [k, n, m]);
  X = permute(X, [1, 3, 2]);
  X = reshape(X, [k * m, n]);
end
