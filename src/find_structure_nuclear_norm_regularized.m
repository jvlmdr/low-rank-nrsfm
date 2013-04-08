% Minimizes
%   1/2 sum_ti ||w_ti - R_t s_ti||^2 + lambda nuclear_norm(S)
%
% structure = find_structure_nuclear_norm_regularized(projections, rotations,
%     lambda, rho, max_iter, mu, tau_incr, tau_decr)
%
% Parameters:
% projections -- 2 x P x F
% structure (optional, or []) -- 3 x P x F. Provided for "warm starting"
% rotations -- 2 x 3 x F
%
% Results:
% structure -- 3 x P x F

function structure = find_structure_nuclear_norm_regularized(projections, ...
    structure, rotations, lambda, rho, tau, rho_max, primal_tol, dual_tol, ...
    max_iter, verbose)
  % Introduce the auxiliary variable X.
  %
  % arg min_{S, X} 1/2 sum_ti ||w_ti - R_t s_ti||^2 + lambda nuclear_norm(X)
  % s.t.  S - X = 0

  if isempty(rho)
    rho = 0;
  end

  P = size(projections, 2);
  F = size(projections, 3);

  if isempty(structure)
    structure = zeros(3, P, F);
  end

  converged = false;
  num_iter = 0;

  % Convex problem, use any initialization.
  S = structure;
  X = structure;
  U = S - X;

  while ~converged && num_iter < max_iter
    prev_S = S;
    prev_X = X;
    prev_U = U;

    % S subproblem. Linear least squares.
    % arg min_{S, X} 1/2 sum_ti ||w_ti - R_t s_ti||^2 + rho/2 ||S - X + U||_F^2
    V = X - U;
    for t = 1:F
      R_t = rotations(:, :, t);
      W_t = projections(:, :, t);
      V_t = V(:, :, t);
      S_t = pinv(R_t' * R_t + rho * eye(3)) * (R_t' * W_t + rho * V_t);
      S(:, :, t) = S_t;
    end

    % X subproblem. Singular value soft thresholding.
    % arg min_{X} nuclear_norm(X) + rho/2 ||S - X + U||_F
    V = S + U;
    V = k_reshape(structure_to_matrix(V), 3);
    if rho == 0
      d = svd(V);
      % lambda / rho == d(1)
      rho = lambda / d(1);
      if verbose
        fprintf('rho (automatic) = %g\n', rho);
      end
    end
    X = singular_value_soft_threshold(V, lambda / rho);
    X = structure_from_matrix(k_unreshape(X, 3));

    % Update multipliers.
    r = S - X;
    U = U + r;
    s = rho * (X - prev_X);

    Y = rho * U;
    norm_r = norm(r(:));
    norm_s = norm(s(:));
    norm_r_rel = norm(r(:)) / max(norm(S(:)), norm(X(:)));
    norm_s_rel = norm(s(:)) / norm(Y(:));
    delta_x = norm(S(:) - prev_S(:)) / norm(prev_S(:));
    delta_z = norm(X(:) - prev_X(:)) / norm(prev_X(:));
    delta_u = norm(U(:) - prev_U(:)) / norm(prev_U(:));
    if verbose
      fprintf(...
          '%6d:  rho:% .2e  r/x:% .2e  s/y:% .2e  dx:% .2e  dz:% .2e  du:% .2e\n', ...
          num_iter, rho, norm_r_rel, norm_s_rel, delta_x, delta_z, delta_u);
    end

    eps_primal = 1e-6 * sqrt(numel(S)) + ...
        primal_tol * max(norm(S(:)), norm(X(:)));
    eps_dual = 1e-6 * sqrt(numel(Y)) + dual_tol * norm(Y(:));

    if norm_r < eps_primal % && norm_s < eps_dual
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

  % Use X to get low rank (to numerical precision).
  structure = X;
end
