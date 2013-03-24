% basis = find_basis(projections, rotations, coeff)
%
% Parameters:
% projections -- 2 x P x F
% rotations -- 2 x 3 x F
% coeff -- K x F
%
% Returns:
% basis -- 3 x K x P

function basis = find_basis(projections, rotations, coeff)
  P = size(projections, 2);
  F = size(projections, 3);
  K = size(coeff, 1);

  % [2, P, F] -> [2, F, P] -> [2F, P]
  W = reshape(permute(projections, [1, 3, 2]), [2 * F, P]);
  % [2, 3, F] -> [2F, 3F]
  R = block_diagonal_cameras(rotations);
  % [K, F] -> [F, K]
  C = coeff';

  B = (R * kron(C, speye(3))) \ W;

  % [3K, P] -> [3, K, P]
  basis = reshape(B, [3, K, P]);
end
