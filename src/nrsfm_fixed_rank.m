% Solves
%   arg min_{R, S} 1/2 sum_t ||W_t - R_t S_t||_F^2
%   s.t. R_t R_t' = I for all t, rank(S) <= K.
%
% [structure, rotations] = nrsfm_fixed_rank(projections, structure, rotations,
%     K, rho1, rho2, max_iter, mu, tau_incr, tau_decr)
%
% Parameters:
% projections -- 2 x P x F
% structure -- 3 x P x F
% rotations -- 2 x 3 x F
% K -- Rank of structure
%
% Returns:
% structure -- 3 x P x F
% rotations -- 2 x 3 x F

function [structure, rotations] = nrsfm_fixed_rank(projections, structure, ...
    rotations, K, rho1, rho2, max_iter, mu, tau_incr, tau_decr)
  % Introduce the auxiliary variables X and Z.
  %
  % arg min_{R, S, X, Z} 1/2 sum_t ||W_t - Z_t X_t||_F^2
  % s.t. R_t R_t' = I for all t, rank(S) <= K.

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

    % S subproblem. Projection on to rank manifold.
    V = X - U1;
    V = structure_to_matrix(V);
    V = k_reshape(V, 3);
    S = project_rank(V, K);
    S = k_unreshape(S, 3);
    S = structure_from_matrix(S);

    % X subproblem. Least squares.
    V = S + U1;
    for t = 1:F
      Z_t = Z(:, :, t);
      W_t = projections(:, :, t);
      V_t = V(:, :, t);
      X_t = (Z_t' * Z_t + rho1 * eye(3)) \ (Z_t' * W_t + rho1 * V_t);
      X(:, :, t) = X_t;
    end

    % R subproblem. Procrustes (nearest orthonormal matrix).
    V = Z - U2;
    for t = 1:F
      R_t = procrustes(eye(2, 3), V(:, :, t));
      R(:, :, t) = R_t(1:2, :);
    end

    % Z subproblem. Least squares.
    for t = 1:F
      W_t = projections(:, :, t);
      Z(:, :, t) = W_t / X(:, :, t);
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

  % Take structure which is low rank.
  structure = X;
  % Take rotations which are on SO(3).
  rotations = R;
end
