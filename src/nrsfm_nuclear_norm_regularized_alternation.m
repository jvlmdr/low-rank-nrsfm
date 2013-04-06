% Solves
%   arg min_{R, X} 1/2 sum_ti ||w_ti - R_t s_ti||^2 + lambda nuclear_norm(X)
%   s.t.  R_t R_t' = I
%
% [structure, rotations] = nrsfm_nuclear_norm_regularizer(projections,
%     structure, rotations, lambda, rho, tau, rho_max, primal_tol, dual_tol,
%     max_iter)
%
% Parameters:
% projections -- 2 x P x F
% structure -- 3 x P x F
% rotations -- 2 x 3 x F
% lambda
%
% Returns:
% structure -- 3 x P x F
% rotations -- 2 x 3 x F

function [structure, rotations] = nrsfm_nuclear_norm_regularized_alternation(...
    projections, structure, rotations, lambda, outer_max_iter, rho, tau, ...
    rho_max, primal_tol, dual_tol, max_iter)
  % Introduce the auxiliary variable S.
  % arg min_{R, X, S} 1/2 sum_ti ||w_ti - R_t x_ti||^2 + lambda nuclear_norm(X)
  % s.t.  R_t R_t' = I,  X - S = 0

  P = size(projections, 2);
  F = size(projections, 3);

  converged = false;
  num_iter = 0;

  while ~converged && num_iter < outer_max_iter
    structure = find_structure_nuclear_norm_regularized(projections, ...
        structure, rotations, lambda, rho, tau, rho_max, primal_tol, ...
        dual_tol, max_iter);

    prev_rotations = rotations;
    for t = 1:F
      W_t = projections(:, :, t);
      X_t = structure(:, :, t);
      R_t = procrustes(X_t', W_t')';
      rotations(:, :, t) = R_t(1:2, :);
    end

    % Rotate entire problem such that rotations are similar to previous set.
    % arg min_{U} 1/2 sum_t ||R_t V - Q_t||
    % [2, 3, F] -> [2, F, 3] -> [2F, 3]
    R = reshape(permute(rotations, [1, 3, 2]), [2 * F, 3]);
    Q = reshape(permute(prev_rotations, [1, 3, 2]), [2 * F, 3]);
    V = procrustes(R, Q);
    for t = 1:F
      rotations(:, :, t) = rotations(:, :, t) * V;
      structure(:, :, t) = V' * structure(:, :, t);
    end

    delta = norm(X(:) - prev_X(:)) / norm(prev_X(:));
    fprintf(...
        '%6d:  rho:% .2e  r:% .2e  s:% .2e  r/x:% .2e  s/y:% .2e  dz:% .2e\n', ...
        num_iter, norm_r, norm_s, norm_r_rel, norm_s_rel, delta);

    if ~converged && num_iter < max_iter
      if rho < rho_max && norm_s < eps_dual && norm_r >= eps_primal
        rho = rho * tau;
      end
    end

    num_iter = num_iter + 1;
  end
end
