% Solves
%   arg min_{R, S} nuclear_norm(S)  s.t.  w_ti = R_t s_ti,  R_t R_t' = I
%
% [structure, rotations] = nrsfm_constrained_nuclear_norm(projections,
%     structure, rotations, rho1, rho2, tau, rho_max, tol, max_iter)
%
% Parameters:
% projections -- 2 x P x F
% structure -- 3 x P x F
% rotations -- 2 x 3 x F
%
% Returns:
% structure -- 3 x P x F
% rotations -- 2 x 3 x F

function [structure, rotations] = nrsfm_constrained_nuclear_norm(...
    projections, structure, rotations, rho1, rho2, tau, rho_max, tol, max_iter)
  % Introduce the auxiliary variables X and Z.
  %
  % arg min_{R, S, X, Z} nuclear_norm(S)
  % s.t. S = X, R = Z, W_t = Z_t X_t, R_t R_t' = I for all t.

  P = size(projections, 2);
  F = size(projections, 3);

  converged = false;
  num_iter = 0;

  S = structure;
  X = structure;
  U1 = S - X;

  R = rotations;
  Z = rotations;
  U2 = R - Z;

  while ~converged && num_iter < max_iter
    prev_X = X;
    prev_Z = Z;

    % R subproblem. Procrustes (nearest orthonormal matrix).
    V = Z - U2;
    for t = 1:F
      %R_t = procrustes(eye(2, 3), V(:, :, t));
      R_t = procrustes(eye(3), V(:, :, t)')';
      R(:, :, t) = R_t(1:2, :);
    end

    % Z subproblem. Linear equations.
    for t = 1:F
      W_t = projections(:, :, t);
      Z(:, :, t) = W_t / X(:, :, t);
    end

    % X subproblem. Nearest point on subspace.
    V = S + U1;
    for t = 1:F
      Z_t = Z(:, :, t);
      W_t = projections(:, :, t);
      V_t = V(:, :, t);
      X_t = V_t + pinv(Z_t) * (W_t - Z_t * V_t);
      X(:, :, t) = X_t;
    end

    % S subproblem. Singular value soft thresholding.
    V = X - U1;
    V = structure_to_matrix(V);
    V = k_reshape(V, 3);
    S = singular_value_soft_threshold(V, 1 / rho1);
    S = k_unreshape(S, 3);
    S = structure_from_matrix(S);

    % Update multipliers and residuals.
    r1 = S - X;
    r2 = R - Z;
    U1 = U1 + r1;
    U2 = U2 + r2;
    s1 = rho1 * (X - prev_X);
    s2 = rho2 * (Z - prev_Z);

    norm_r1 = norm(r1(:));
    norm_r2 = norm(r2(:));
    norm_s1 = norm(s1(:));
    norm_s2 = norm(s2(:));

    fprintf('%6d %11.3g %11.3g %11.3g %11.3g %11.3g %11.3g\n', num_iter, ...
        rho1, norm_r1, norm_s1, rho2, norm_r2, norm_s2);

    Y1 = rho1 * U1;
    Y2 = rho2 * U2;

    eps_primal1 = sqrt(numel(S)) * 1e-6 + tol * max(norm(S(:)), norm(X(:)));
    eps_dual1 = sqrt(numel(S)) * 1e-6 + tol * norm(Y1(:));

    eps_primal2 = sqrt(numel(R)) * 1e-6 + tol * max(norm(R(:)), norm(Z(:)));
    eps_dual2 = sqrt(numel(R)) * 1e-6 + tol * norm(Y2(:));

    if norm_r1 < eps_primal1 && norm_r2 < eps_primal2
      converged = true;
    end

    if ~converged && num_iter < max_iter
      if rho1 < rho_max && rho2 < rho_max
        if norm_r1 < norm_r2
          rho2 = rho2 * tau;
        else
          rho1 = rho1 * tau;
        end
      else
        % At most one of { rho1, rho2 } can be increased.
        if rho1 < rho_max
          rho1 = rho1 * tau;
        elseif rho2 < rho_max
          rho2 = rho2 * tau;
        end
      end
    end

    U1 = 1 / rho1 * Y1;
    U2 = 1 / rho2 * Y2;
    num_iter = num_iter + 1;
  end

  structure = S;
  rotations = R;
end
