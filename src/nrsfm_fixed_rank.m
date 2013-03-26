% Solves
%   arg min_{R, S} 1/2 sum_t ||Wt - Rt St||_F^2
%   s.t. Rt Rt' = I for all t, rank(S) <= K.
%
% Parameters:
% W -- 2 x P x F
% R_init -- 2 x 3 x F
% K -- Rank of structure
%
% Returns:
% S -- 3 x P x F
% R -- 2 x 3 x F

function [S, R] = nrsfm_fixed_rank(W, R_init, S_init, K, rho1, rho2, ...
    max_iter, mu, tau_incr, tau_decr)
  % Introduce the auxiliary variables X and Z.
  %
  % arg min_{R, S, X, Z} 1/2 sum_t ||Wt - Zt Xt||_F^2
  % s.t. Rt Rt' = I for all t, rank(S) <= K.

  P = size(W, 2);
  F = size(W, 3);

  converged = false;
  num_iter = 0;

  S = S_init;
  X = S_init;
  U1 = S - X;

  R = R_init;
  Z = R_init;
  U2 = R - Z;

  while ~converged && num_iter < max_iter
    prev_X = X;
    prev_Z = Z;

    % S subproblem. Projection on to rank manifold.
    V = X - U1;
    V = reshape(V, [3 * P, F]);
    S = project_rank(V, K);
    S = reshape(S, [3, P, F]);

    % X subproblem. Least squares.
    V = S + U1;
    for t = 1:F
      Zt = Z(:, :, t);
      Wt = W(:, :, t);
      Vt = V(:, :, t);
      Xt = (Zt' * Zt + rho1 * eye(3)) \ (Zt' * Wt + rho1 * Vt);
      X(:, :, t) = Xt;
    end

    % R subproblem. Procrustes (nearest orthonormal matrix).
    V = Z - U2;
    for t = 1:F
      Rt = procrustes(eye(2, 3), V(:, :, t));
      R(:, :, t) = Rt(1:2, :);
    end

    % Z subproblem. Least squares.
    for t = 1:F
      Z(:, :, t) = W(:, :, t) / X(:, :, t);
    end

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

    fprintf('%12d %12g %12g %12g %12g %12g %12g\n', num_iter, rho1, norm_r1, ...
        norm_s1, rho2, norm_r2, norm_s2);

    % Only consider changing rho in first half of iterations.
    if num_iter > 0 && num_iter < max_iter / 2
      if norm_r1 ~= 0 && norm_s1 ~= 0
        if norm_r1 > mu * norm_s1
          rho1 = rho1 * tau_incr;
        elseif norm_s1 > mu * norm_r1
          rho1 = rho1 / tau_decr;
        end
      end

      if norm_r2 ~= 0 && norm_s2 ~= 0
        if norm_r2 > mu * norm_s2
          rho2 = rho2 * tau_incr;
        elseif norm_s2 > mu * norm_r2
          rho2 = rho2 / tau_decr;
        end
      end
    end

    num_iter = num_iter + 1;
  end
end
