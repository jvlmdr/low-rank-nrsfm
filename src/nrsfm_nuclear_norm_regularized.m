% Solves
%   arg min_{R, X} 1/2 sum_ti ||w_ti - R_t s_ti||^2 + lambda nuclear_norm(X)
%   s.t.  R_t R_t' = I
%
% [structure, rotations] = nrsfm_nuclear_norm_regularized(projections,
%     structure, rotations, lambda, rho1, rho2, tau, rho_max, primal_tol,
%     dual_tol, max_iter)
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

function [structure, rotations] = nrsfm_nuclear_norm_regularized(...
    projections, structure, rotations, lambda, rho1, rho2, tau, rho_max, ...
    primal_tol, dual_tol, max_iter)
  % Introduce the auxiliary variable S.
  % arg min_{R, X, S} 1/2 sum_ti ||w_ti - R_t x_ti||^2 + lambda nuclear_norm(X)
  % s.t.  R_t R_t' = I,  X - S = 0

  P = size(projections, 2);
  F = size(projections, 3);

  f_init = 1/2 * projection_error(projections, structure, rotations) ^ 2 + ...
      lambda * nuclear_norm(k_reshape(structure_to_matrix(structure), 3));

  converged = false;
  num_iter = 0;

  R = rotations;
  Z = rotations;
  X = structure;
  S = structure;
  U1 = S - X;
  U2 = R - Z;

  while ~converged && num_iter < max_iter
    prev_X = X;
    prev_Z = Z;

    % S subproblem. Linear least squares.
    % arg min_{S} 1/2 sum_ti ||w_ti - R_t s_ti||^2 + rho1/2 ||S - X + U1||_F^2
    V = X - U1;
    for t = 1:F
      R_t = R(:, :, t);
      W_t = projections(:, :, t);
      V_t = V(:, :, t);
      X_t = (R_t' * R_t + rho1 * eye(3)) \ (R_t' * W_t + rho1 * V_t);
      S(:, :, t) = X_t;
    end

    % X subproblem. Singular value soft thresholding.
    % arg min_{X} lambda nuclear_norm(X) + rho1/2 ||S - X + U1||_F^2
    V = S + U1;
    V = structure_to_matrix(V);
    V = k_reshape(V, 3);
    if rho1 == 0
      d = svd(V);
      rho1 = lambda / (1e-2 * d(1));
    end
    X = singular_value_soft_threshold(V, lambda / rho1);
    X = k_unreshape(X, 3);
    X = structure_from_matrix(X);

    % R subproblem. Procrustes.
    % arg min_{R} ||R - Z + U2||_F^2  s.t.  R_t R_t' = I
    prev_R = R;
    V = Z - U2;
    for t = 1:F
      R_t = procrustes(eye(2, 3), Z(:, :, t));
      R(:, :, t) = R_t(1:2, :);
    end

    % Z subproblem. Linear least squares.
    % arg min_{Z} 1/2 ||W - Z X||_F^2 + rho2/2 ||R - Z + U2||_F^2
    V = R + U2;
    for t = 1:F
      W_t = projections(:, :, t);
      X_t = X(:, :, t);
      R_t = R(:, :, t);
      V_t = V(:, :, t);
      Z_t = (X_t * X_t' + rho2 * eye(3)) \ (X_t * W_t' + rho2 * V_t');
      Z(:, :, t) = Z_t';
    end

    % Rotate entire problem so that R changes minimally.
    % arg min_{V} 1/2 sum_t ||R_t V - Q_t||
    % [2, 3, F] -> [2, F, 3] -> [2F, 3]
    A = reshape(permute(R, [1, 3, 2]), [2 * F, 3]);
    B = reshape(permute(prev_R, [1, 3, 2]), [2 * F, 3]);
    V = procrustes(A, B);
    for t = 1:F
      S(:, :, t) = V' * S(:, :, t);
      X(:, :, t) = V' * X(:, :, t);
      R(:, :, t) = R(:, :, t) * V;
      Z(:, :, t) = Z(:, :, t) * V;
    end

    % Update multipliers and residuals.
    r1 = S - X;
    s1 = rho1 * (X - prev_X);
    U1 = U1 + r1;
    Y1 = rho1 * U1;

    r2 = R - Z;
    s2 = rho2 * (Z - prev_Z);
    U2 = U2 + r2;
    Y2 = rho2 * U2;

    norm_r1 = norm(r1(:));
    norm_s1 = norm(s1(:));
    norm_r1_rel = norm(r1(:)) / max(norm(S(:)), norm(X(:)));
    norm_s1_rel = norm(s1(:)) / norm(Y1(:));

    norm_r2 = norm(r2(:));
    norm_s2 = norm(s2(:));
    norm_r2_rel = norm(r2(:)) / max(norm(R(:)), norm(Z(:)));
    norm_s2_rel = norm(s2(:)) / norm(Y2(:));

    fprintf(...
        '%6d:  rho1:% .2e  r/x:% .2e  s/y:% .2e  rho2:% .2e  r/x:% .2e  s/y:% .2e\n', ...
        num_iter, rho1, norm_r1_rel, norm_s1_rel, rho2, norm_r2_rel, ...
        norm_s2_rel);

    primal_eps1 = 1e-6 + primal_tol * max(norm(S(:)), norm(X(:)));
    dual_eps1 = 1e-6 + dual_tol * norm(Y1(:));
    primal_eps2 = 1e-6 + primal_tol * max(norm(R(:)), norm(Z(:)));
    dual_eps2 = 1e-6 + dual_tol * norm(Y2(:));

    if norm_r1 <= primal_eps1 && norm_r2 <= primal_eps2
      converged = true;
    end

    if ~converged && num_iter < max_iter
      if norm_s1 < dual_eps1 && norm_s2 < dual_eps2
        if rho1 < rho_max && norm_r1 > primal_eps1
          rho1 = rho1 * tau;
        end
        if rho2 < rho_max && norm_r2 > primal_eps2
          rho2 = rho2 * tau;
        end
      end
    end

    U1 = 1 / rho1 * Y1;
    U2 = 1 / rho2 * Y2;
    num_iter = num_iter + 1;
  end

  structure = X;
  rotations = R;

  f = projection_error(projections, structure, rotations) ^ 2 + ...
      lambda * nuclear_norm(k_reshape(structure_to_matrix(structure), 3));

  fprintf('Initial objective: %g\n', f_init);
  fprintf('Final objective: %g\n', f);
end
