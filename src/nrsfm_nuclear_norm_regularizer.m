% Solves
% arg min_{R, S} 1/2 sum_ti ||w_ti - R_t s_ti||^2 + lambda nuclear_norm(S)
% s.t.  R_t R_t' = I
%
% [structure, rotations] = nrsfm_nuclear_norm_regularizer(projections,
%     rotations, basis, rho, max_iter, mu, tau_incr, tau_decr)
%
% Parameters:
% projections -- 2 x P x F
% rotations -- 2 x 3 x F
% lambda
%
% Returns:
% structure -- 3 x P x F
% rotations -- 2 x 3 x F

function [structure, rotations] = nrsfm_nuclear_norm_regularizer(...
    projections, rotations, lambda, rho, max_iter, mu, tau_incr, tau_decr)
  % Introduce the auxiliary variable X.
  % arg min_{R, S, X} 1/2 sum_ti ||w_ti - R_t x_ti||^2 + lambda nuclear_norm(S)
  % s.t.  R_t R_t' = I,  S - X = 0

  P = size(projections, 2);
  F = size(projections, 3);

  converged = false;
  num_iter = 0;

  R = rotations;
  S = zeros(3, P, F);
  X = zeros(3, P, F);
  U = S - X;

  while ~converged && num_iter < max_iter
    prev_X = X;

    % X subproblem. Linear least squares.
    % arg min_{X} 1/2 sum_ti ||w_ti - R_t x_ti||^2 + rho/2 ||S - X + U||_F^2
    V = S + U;
    for t = 1:F
      R_t = R(:, :, t);
      W_t = projections(:, :, t);
      V_t = V(:, :, t);
      X_t = (R_t' * R_t + rho * eye(3)) \ (R_t' * W_t + rho * V_t);
      X(:, :, t) = X_t;
    end

    % S subproblem. Singular value soft thresholding.
    % arg min_{S} lambda nuclear_norm(S) + rho/2 ||S - X + U||_F^2
    V = X - U;
    V = reshape(V, [3 * P, F]);
    S = singular_value_soft_threshold(V, lambda / rho);
    S = reshape(S, [3, P, F]);

    % R subproblem. Orthogonal Procrustes problem
    % arg min_{R} 1/2 sum_t ||W_t - R_t X_t||^2  s.t.  R_t R_t' = I
    for t = 1:F
      W_t = projections(:, :, t);
      X_t = X(:, :, t);
      R_t = procrustes(X_t', W_t')';
      R(:, :, t) = R_t(1:2, :);
    end

    % Update multipliers and residuals.
    r = S - X;
    U = U + r;
    s = rho * (X - prev_X);

    norm_r = norm(r(:));
    norm_s = norm(s(:));

    fprintf('%6d: %12g %12g %12g\n', num_iter, rho, norm_r, norm_s);

    % Only consider changing rho in first half of iterations.
    if num_iter > 0 && num_iter < max_iter / 2
      if norm_r ~= 0 && norm_s ~= 0
        if norm_r > mu * norm_s
          rho = rho * tau_incr;
        elseif norm_s > mu * norm_r
          rho = rho / tau_decr;
        end
      end
    end

    num_iter = num_iter + 1;
  end

  structure = S;
  rotations = R;
end
