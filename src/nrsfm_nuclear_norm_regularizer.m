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

function [structure, rotations] = nrsfm_nuclear_norm_regularizer(...
    projections, structure, rotations, lambda, rho, tau, rho_max, ...
    primal_tol, dual_tol, max_iter)
  % Introduce the auxiliary variable S.
  % arg min_{R, X, S} 1/2 sum_ti ||w_ti - R_t x_ti||^2 + lambda nuclear_norm(X)
  % s.t.  R_t R_t' = I,  X - S = 0

  P = size(projections, 2);
  F = size(projections, 3);

  if isempty(rho)
    rho = 0;
  end

  converged = false;
  num_iter = 0;

  R = rotations;
  X = structure;
  S = structure;
  U = S - X;

  while ~converged && num_iter < max_iter
    prev_X = X;

    % X subproblem. Singular value soft thresholding.
    % arg min_{X} lambda nuclear_norm(X) + rho/2 ||S - X + U||_F^2
    V = S + U;
    V = structure_to_matrix(V);
    V = k_reshape(V, 3);
    if rho == 0
      d = svd(V);
      rho = lambda / (1e-2 * d(1));
    end
    X = singular_value_soft_threshold(V, lambda / rho);
    X = k_unreshape(X, 3);
    X = structure_from_matrix(X);

    % S subproblem. Linear least squares.
    % arg min_{S} 1/2 sum_ti ||w_ti - R_t s_ti||^2 + rho/2 ||S - X + U||_F^2
    V = X - U;
    for t = 1:F
      R_t = R(:, :, t);
      W_t = projections(:, :, t);
      V_t = V(:, :, t);
      X_t = (R_t' * R_t + rho * eye(3)) \ (R_t' * W_t + rho * V_t);
      S(:, :, t) = X_t;
    end

    % R subproblem. Orthogonal Procrustes problem
    % arg min_{R} 1/2 sum_t ||W_t - R_t X_t||^2  s.t.  R_t R_t' = I
    Q = R;
    for t = 1:F
      W_t = projections(:, :, t);
      X_t = S(:, :, t);
      R_t = procrustes(X_t', W_t')';
      R(:, :, t) = R_t(1:2, :);
    end

    % Rotate entire problem such that rotations are similar to previous set.
    % arg min_{U} 1/2 sum_t ||R_t V - Q_t||
    % [2, 3, F] -> [2, F, 3] -> [2F, 3]
    RR = reshape(permute(R, [1, 3, 2]), [2 * F, 3]);
    QQ = reshape(permute(Q, [1, 3, 2]), [2 * F, 3]);
    V = procrustes(RR, QQ);
    for t = 1:F
      R(:, :, t) = R(:, :, t) * V;
      S(:, :, t) = V' * S(:, :, t);
      X(:, :, t) = V' * X(:, :, t);
    end

    % Update multipliers and residuals.
    r = S - X;
    U = U + r;
    s = X - prev_X;

    Y = rho * U;
    eps_primal = sqrt(numel(S)) * 1e-6 + ...
        primal_tol * max(norm(S(:)), norm(X(:)));
    eps_dual = sqrt(numel(U)) * 1e-6 + dual_tol * norm(U(:));

    norm_r = norm(r(:));
    norm_s = norm(s(:));
    norm_r_rel = norm(r(:)) / max(norm(S(:)), norm(X(:)));
    norm_s_rel = norm(s(:)) / norm(U(:));
    delta = norm(X(:) - prev_X(:)) / max(norm(X(:)), norm(prev_X(:)));
    fprintf('%6d:  rho:% .2e  r/x:% .2e  s/y:% .2e  dz:% .2e\n', num_iter, ...
        rho, norm_r_rel, norm_s_rel, delta);

    if norm_r < eps_primal
      converged = true;
    end

    if ~converged && num_iter < max_iter
      if rho < rho_max && norm_s < eps_dual && norm_r >= eps_primal
        rho = rho * tau;
      end
    end

    U = 1 / rho * Y;
    num_iter = num_iter + 1;
  end

  structure = X;
  rotations = R;
end
