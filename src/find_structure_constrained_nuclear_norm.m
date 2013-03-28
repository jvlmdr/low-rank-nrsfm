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
%
% Results:
% structure -- 3 x P x F

function structure = find_structure_constrained_nuclear_norm(projections, ...
    rotations, rho, max_iter, mu, tau_incr, tau_decr)
  % Introduce the auxiliary variable Z.
  %
  % arg min_{S, Z} nuclear_norm(Z)
  % s.t.  S = Z,  w_ti = R_t z_ti

  P = size(projections, 2);
  F = size(projections, 3);

  converged = false;
  num_iter = 0;

  % Convex problem, use any initialization.
  S = zeros(3, P, F);
  Z = zeros(3, P, F);
  U = S - Z;

  while ~converged && num_iter < max_iter
    prev_Z = Z;

    % S subproblem. Nearest point on subspace.
    % arg min_{S} ||S - Z + U||_F  s.t.  W_t = R_t S_t
    V = Z - U;
    for t = 1:F
      R_t = rotations(:, :, t);
      W_t = projections(:, :, t);
      V_t = V(:, :, t);
      % For some reason this is very wrong:
      % S_t = V_t + R_t \ (W_t - R_t * V_t);
      S_t = V_t + pinv(R_t) * (W_t - R_t * V_t);
      S(:, :, t) = S_t;
    end

    % Z subproblem. Singular value soft thresholding.
    % arg min_{Z} nuclear_norm(Z) + rho/2 ||S - Z + U||_F
    V = S + U;
    V = structure_to_matrix(V);
    V = k_reshape(V, 3);
    Z = singular_value_soft_threshold(V, 1 / rho);
    Z = k_unreshape(Z, 3);
    Z = structure_from_matrix(Z);

    % Update multipliers.
    r = S - Z;
    U = U + r;
    s = rho * (Z - prev_Z);

    norm_r = norm(r(:));
    norm_s = norm(s(:));

    fprintf('%6d: %8g %8g %8g\n', num_iter, rho, norm_r, norm_s);

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

  % Use Z to get low rank (to numerical precision).
  structure = Z;
end
