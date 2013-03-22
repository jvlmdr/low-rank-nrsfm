% Parameters:
% points -- 3 x P x F
% K -- Basis size
%
% Returns:
% basis -- 3 x K x P
% coeff -- K x F

function [basis, coeff] = factorize_structure(points, K)
  % [3, P, F] -> [3P, F]
  S = points;
  S = reshape(S, [3 * P, F]);

  [U, S, V] = svd(S, 'econ');
  B = U(:, 1:K);
  C = V(:, 1:K) * D(1:K, 1:K);

  % [3P, K] -> [3, P, K] -> [3, K, P]
  basis = B;
  basis = reshape(basis, [3, P, K]);
  basis = permute(basis, [1, 3, 2]);
  % [F, K] -> [K, F]
  coeff = C';
end
