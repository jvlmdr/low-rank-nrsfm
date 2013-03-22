% [structure, basis, coeff] =
%   find_structure_nullspace(projections, rotations, K)
%
% Parameters:
% projections -- 2 x P x F
% rotations -- 2 x 3 x F
% K -- Shape basis size
%
% Returns:
% structure -- 3 x P x F
% basis -- 3 x K x P
% coeff -- K x F

function [structure, basis, coeff] = find_structure_nullspace(projections, ...
    rotations, K)
  P = size(projections, 2);
  F = size(projections, 3);

  [M_hat, B_hat, W] = factorize_projections(projections, K);

  R = block_diagonal_cameras(rotations);
  [G, C] = find_corrective_matrix_nullspace(M_hat, R);
  B = inv(G) * B_hat;
  S = kron(C, eye(3)) * B;

  % [3F, P] -> [3, F, P] -> [3, P, F]
  structure = S;
  structure = reshape(structure, [3, F, P]);
  structure = permute(structure, [1, 3, 2]);
  % [3K, P] -> [3, K, P]
  basis = B;
  basis = reshape(basis, [3, K, P]);
  % [F, K] -> [K, F]
  coeff = C';
end
