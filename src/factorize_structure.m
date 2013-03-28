% Parameters:
% structure -- 3 x P x F
% K -- Basis size
%
% Returns:
% basis -- 3 x K x P
% coeff -- K x F

function [basis, coeff] = factorize_structure(structure, K)
  F = size(structure, 3);
  P = size(structure, 2);

  S = k_reshape(structure_to_matrix(structure), 3);

  [U, D, V] = svd(S, 'econ');
  B = U(:, 1:K);
  C = V(:, 1:K) * D(1:K, 1:K);

  B = k_unreshape(B, 3);
  basis = basis_from_matrix(B);
  coeff = C';
end
