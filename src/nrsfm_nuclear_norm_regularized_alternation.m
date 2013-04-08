% Solves
%   arg min_{R, S} 1/2 sum_ti ||w_ti - R_t s_ti||^2 + lambda nuclear_norm(S)
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
    projections, structure, rotations, lambda, outer_tol, outer_max_iter, ...
    rho, tau, rho_max, primal_tol, dual_tol, max_iter)
  P = size(projections, 2);
  F = size(projections, 3);

  converged = false;
  num_iter = 0;
  y = 1/2 * projection_error(projections, structure, rotations) ^ 2 + ...
      lambda * nuclear_norm(k_reshape(structure_to_matrix(structure), 3));
  y_init = y;

  while ~converged && num_iter < outer_max_iter
    prev_structure = structure;
    prev_rotations = rotations;
    prev_y = y;

    structure = find_structure_nuclear_norm_regularized(projections, ...
        structure, rotations, lambda, rho, tau, rho_max, primal_tol, ...
        dual_tol, max_iter, false);
    y1 = 1/2 * projection_error(projections, structure, rotations) ^ 2 + ...
        lambda * nuclear_norm(k_reshape(structure_to_matrix(structure), 3));

    rotations = find_cameras(projections, structure);

    %for t = 1:F
    %  W_t = projections(:, :, t);
    %  X_t = structure(:, :, t);
    %  R_t = procrustes(X_t', W_t')';
    %  keyboard;
    %  rotations(:, :, t) = R_t(1:2, :);
    %end
    y2 = 1/2 * projection_error(projections, structure, rotations) ^ 2 + ...
        lambda * nuclear_norm(k_reshape(structure_to_matrix(structure), 3));

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
    y3 = 1/2 * projection_error(projections, structure, rotations) ^ 2 + ...
        lambda * nuclear_norm(k_reshape(structure_to_matrix(structure), 3));
    assert(abs(y3 - y2) < 1e-6 + 1e-6 * abs(max(y2, y3)));
    y = y2;

    dy = abs(prev_y - y);
    dx = norm(structure(:) - prev_structure(:));
    dy_threshold = 1e-6 + outer_tol * prev_y;
    dx_threshold = 1e-6 + outer_tol * norm(prev_structure(:));

    dy_rel = dy / y;
    dx_rel = dx / norm(prev_structure(:));
    fprintf(...
        '%6d:  y1:% .4e y2:% .4e dy:% .2e  dx:% .2e  dy/y:% .2e  dx/x:% .2e\n', ...
        num_iter, y1, y2, dy, dx, dy_rel, dx_rel);

    if dy < dy_threshold
      converged = true;
    elseif dx < dx_threshold
      converged = true;
    end

    num_iter = num_iter + 1;
  end

  fprintf('Initial objective: %g\n', y_init);
  fprintf('Final objective: %g\n', y);
end
