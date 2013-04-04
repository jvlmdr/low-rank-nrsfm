% Minimizes
%   sum_ti ||w_ti - R_t s_ti||^2  s.t.  rank(S) <= K
% by solving
%   1/2 sum_ti ||w_ti - R_t s_ti||^2 + lambda nuclear_norm(S)
% for varying lambda.
%
% [structure, basis, coeff] = find_structure_nuclear_norm_sweep(projections,
%     rotations, K, rho, max_iter, mu, tau_incr, tau_decr)
%
% Parameters:
% projections -- 2 x P x F
% rotations -- 2 x 3 x F
%
% Results:
% structure -- 3 x P x F
% basis -- 3 x K x P
% coeff -- K x F

function [structure, basis, coeff] = find_structure_nuclear_norm_sweep(...
    projections, rotations, K, rho, max_iter, tol, tau, rho_max)
  F = size(projections, 3);
  P = size(projections, 2);

  lambdas = logspace(-6, 6, 13);
  num_lambdas = length(lambdas);

  structures = zeros(3, P, F, num_lambdas);
  bases = zeros(3, K, P, num_lambdas);
  coeffs = zeros(K, F, num_lambdas);

  for i = 1:num_lambdas
    structure = find_structure_nuclear_norm_regularized(projections, ...
        rotations, lambdas(i), rho, max_iter, tol, tau, rho_max);
    [basis, coeff] = factorize_structure(structure, K);
    structure = compose_structure(basis, coeff);

    bases(:, :, :, i) = basis;
    coeffs(:, :, i) = coeff;
    structures(:, :, :, i) = structure;
  end

  W = projections_to_matrix(projections);
  R = block_diagonal_cameras(rotations);

  residuals = zeros(num_lambdas, 1);
  for i = 1:num_lambdas
    S = structure_to_matrix(structures(:, :, :, i));
    residuals(i) = norm(W - R * S, 'fro');
  end

  [residual, i] = min(residuals);
  structure = structures(:, :, :, i);
  basis = bases(:, :, :, i);
  coeff = coeffs(:, :, i);
end
