% Minimizes
%   nuclear_norm(S)
% subject to
%   w_ti = R_t s_ti.
%
% structure = find_structure_constrained_nuclear_norm(projections, rotations,
%     rho, max_iter, mu, tau_incr, tau_decr)
%
% Parameters:
% projections -- 2 x P x F
% rotations -- 2 x 3 x F
% rho -- Can be []
%
% Results:
% structure -- 3 x P x F

function structure = find_structure_constrained_nuclear_norm(projections, ...
    rotations, rho, tau, rho_max, primal_tol, dual_tol, max_iter, verbose)
  % Introduce the auxiliary variable X.
  %
  % arg min_{S, X} nuclear_norm(X)
  % s.t.  S - X = 0,  w_ti = R_t z_ti

  P = size(projections, 2);
  F = size(projections, 3);

  converged = false;
  num_iter = 0;

  % Convex problem, use any initialization.
  S = zeros(3, P, F);
  X = zeros(3, P, F);
  U = S - X;

  while ~converged && num_iter < max_iter
    prev_X = X;

    % S subproblem. Nearest point on subspace.
    % arg min_{S} ||S - X + U||_F  s.t.  W_t = R_t S_t
    V = X - U;
    for t = 1:F
      R_t = rotations(:, :, t);
      W_t = projections(:, :, t);
      V_t = V(:, :, t);
      % For some reason this is very wrong:
      % S_t = V_t + R_t \ (W_t - R_t * V_t);
      S_t = V_t + pinv(R_t) * (W_t - R_t * V_t);
      S(:, :, t) = S_t;
    end

    % X subproblem. Singular value soft thresholding.
    % arg min_{X} nuclear_norm(X) + rho/2 ||S - X + U||_F
    V = S + U;
    V = structure_to_matrix(V);
    V = k_reshape(V, 3);
    % Set initial rho automatically if empty.
    if isempty(rho)
      d = svd(V);
      rho = 1 / d(1);
      fprintf('rho (automatic) = %g\n', rho);
    end
    X = singular_value_soft_threshold(V, 1 / rho);
    X = k_unreshape(X, 3);
    X = structure_from_matrix(X);

    % Update multipliers.
    r = S - X;
    U = U + r;
    s = rho * (X - prev_X);

    Y = rho * U;
    norm_r = norm(r(:));
    norm_s = norm(s(:));
    norm_r_rel = norm(r(:)) / max(norm(S(:)), norm(X(:)));
    norm_s_rel = norm(s(:)) / norm(Y(:));
    if verbose
      fprintf('%6d:  rho:% .2e  r/x:% .2e  s/y:% .2e\n', num_iter, rho, ...
          norm_r_rel, norm_s_rel);
    end

    eps_primal = 1e-6 + primal_tol * max(norm(S(:)), norm(X(:)));
    eps_dual = 1e-6 + dual_tol * norm(Y(:));

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
