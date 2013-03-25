% Solves NRSFM using method of Xiao and Kanade (2004).
%
% [structure, rotations, basis, coeff] = nrsfm_basis_constraints(projections, K)
%
% Parameters:
% projections -- 2 x P x F
% K -- Basis size
%
% Returns:
% structure -- 3 x P x F
% rotations -- 2 x 3 x F
% basis -- 3 x K x P
% coeff -- K x F
%
% Projections do not need to have centroid removed already.

function [structure, rotations, basis, coeff] = nrsfm_basis_constraints(...
    projections, K)
  F = size(projections, 3);

  [M_hat, B_hat, W] = factorize_projections(projections, K);

  % Use a random subset of K frames for basis.
  subset = randperm(F);
  subset = subset(1:K);

  % Solve for corrective matrix.
  [G, rotations, coeff] = find_corrective_matrix_basis_constraints(M_hat, subset);
  % Recover structure.
  B = inv(G) * B_hat;

  basis = basis_from_matrix(B);
  structure = compose_structure(basis, coeff);
end
