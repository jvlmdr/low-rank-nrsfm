% [structure, rotations, basis, coeff] =
%   nrsfm_nullspace_alternation(projections, rotations, basis, num_iter)
%
% Parameters:
% projections -- 2 x P x F
% rotations -- 2 x 3 x F
% basis -- 3 x K x P
% num_iter -- Number of iterations of alternation to do
%
% Returns:
% structure -- 3 x P x F
% rotations -- 2 x 3 x F
% basis -- 3 x K x P
% coeff -- K x F

function [structure, rotations, basis, coeff] = nrsfm_alternation(...
    projections, rotations, basis, num_iter)
  for i = 1:num_iter
    coeff = find_coefficients(projections, rotations, basis);
    err1 = objective(projections, rotations, basis, coeff);

    structure = compose_structure(basis, coeff);
    rotations = find_non_worse_cameras(projections, structure, rotations);
    err2 = objective(projections, rotations, basis, coeff);

    basis = find_basis(projections, rotations, coeff);
    err3 = objective(projections, rotations, basis, coeff);

    fprintf('%6d: %8g %8g %8g\n', i, err1, err2, err3);
  end
end

function err = objective(projections, rotations, basis, coeff)
  P = size(projections, 2);
  F = size(projections, 3);

  structure = compose_structure(basis, coeff);

  % [2, P, F] -> [2, F, P] -> [2F, P]
  W = reshape(permute(projections, [1, 3, 2]), [2 * F, P]);
  % [2, 3, F] -> [2F, 3F]
  R = block_diagonal_cameras(rotations);
  % [3, P, F] -> [3, F, P] -> [3F, P]
  S = reshape(permute(structure, [1, 3, 2]), [3 * F, P]);

  err = norm(W - R * S, 'fro');
end
