% Solves
%   arg min_{R, S} nuclear_norm(S)  s.t.  w_ti = R_t s_ti,  R_t R_t' = I
%
% [structure, rotations] = nrsfm_constrained_nuclear_norm_biaffine(projections,
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

function [structure, rotations] = nrsfm_constrained_nuclear_norm_biaffine(...
    projections, structure, rotations, rho1, rho2, tau1, tau2, rho_max1, ...
    rho_max2, tol, max_iter)
  % Introduce the auxiliary variables X and Z.
  %
  % arg min_{R, S, X, Z} nuclear_norm(S)
  % s.t. S = X, R = Z, W_t = Z_t X_t, R_t R_t' = I for all t.

  P = size(projections, 2);
  F = size(projections, 3);

  converged = false;
  num_iter = 0;

  W = projections;
  R = rotations;
  S = structure;
  X = structure;

  RR = block_diagonal_cameras(R);
  XX = structure_to_matrix(X);
  U1 = projections_from_matrix(RR * XX) - W;
  U2 = S - X;

  while ~converged && num_iter < max_iter
    prev_S = S;
    prev_X = X;

    % R subproblem. Procrustes (nearest orthonormal matrix).
    % min_R  sum_t ||R_t X_t - W_t + U_t||_F  s.t.  R_t R_t' = I
    V = W - U1;
    for t = 1:F
      R_t = procrustes(X(:, :, t)', V(:, :, t)')';
      R(:, :, t) = R_t(1:2, :);
    end

    % S subproblem. Singular value soft thresholding.
    % min_S nuclear_norm(S) + rho/2 ||S - X + U||_F^2
    V = X - U2;
    V = structure_to_matrix(V);
    V = k_reshape(V, 3);
    S = singular_value_soft_threshold(V, 1 / rho1);
    S = k_unreshape(S, 3);
    S = structure_from_matrix(S);

    % X subproblem. Linear least squares.
    % min_X rho1/2 ||R S - W + U1||_F^2 + rho2/2 ||S - X + U2||_F^2 
    % min_X rho1/2 ||R S - V1||_F^2 + rho2/2 ||V2 - X||_F^2 
    V1 = W - U1;
    V2 = S + U2;
    for t = 1:F
      X_t = (rho1 * R_t' * R_t + rho2 * eye(3)) \ ...
          (rho1 * R_t' * V1(:, :, t) + rho2 * V2(:, :, t));
      X(:, :, t) = X_t;
    end

    % Update multipliers and residuals.
    RR = block_diagonal_cameras(R);
    XX = structure_to_matrix(X);
    r1 = projections_from_matrix(RR * XX) - W;
    r2 = S - X;

    U1 = U1 + r1;
    U2 = U2 + r2;
    s1 = rho1 * (S - prev_S);
    s2 = rho2 * (X - prev_X);

    norm_r1 = norm(r1(:));
    norm_r2 = norm(r2(:));
    norm_s1 = norm(s1(:));
    norm_s2 = norm(s2(:));

    fprintf('%6d %11.3g %11.3g %11.3g %11.3g %11.3g %11.3g\n', num_iter, ...
        rho1, norm_r1, norm_s1, rho2, norm_r2, norm_s2);

    Y1 = rho1 * U1;
    Y2 = rho2 * U2;

    % Residual for g(R, S) = W.
    eps_primal1 = 1e-6 * sqrt(max(numel(R), numel(S))) + ...
        tol * max(max(norm(R(:)), norm(S(:))), norm(W(:)));
    eps_dual1 = 1e-6 * sqrt(numel(Y1)) + tol * norm(Y1(:));
    % Residual for S - X = 0.
    eps_primal2 = 1e-6 * sqrt(numel(S)) + tol * max(norm(S(:)), norm(X(:)));
    eps_dual2 = 1e-6 * sqrt(numel(Y2)) + tol * norm(Y2(:));

    if norm_r1 < eps_primal1 && norm_r2 < eps_primal2
      converged = true;
    end

    if ~converged && num_iter < max_iter
      if norm_r1 >= eps_primal1 && rho1 < rho_max1
        rho1 = rho1 * tau1;
      end

      if norm_r2 >= eps_primal2 && rho2 < rho_max2
        rho2 = rho2 * tau2;
      end
    end

    U1 = 1 / rho1 * Y1;
    U2 = 1 / rho2 * Y2;
    num_iter = num_iter + 1;
  end

  structure = S;
  rotations = R;
end
