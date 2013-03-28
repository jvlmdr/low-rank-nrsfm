% Minimizes
%   1/2 sum_ti ||w_ti - R_t s_ti||^2 + lambda nuclear_norm(S)
%
% structure = find_structure_nuclear_norm_regularized(projections, rotations,
%     lambda, rho, max_iter, mu, tau_incr, tau_decr)
%
% Parameters:
% projections -- 2 x P x F
% rotations -- 2 x 3 x F
%
% Results:
% structure -- 3 x P x F

function structure = find_structure_nuclear_norm_regularized(projections, ...
    rotations, lambda, rho, max_iter, mu, tau_incr, tau_decr)
  % Introduce the auxiliary variable Z.
  %
  % arg min_{S, Z} 1/2 sum_ti ||w_ti - R_t s_ti||^2 + lambda nuclear_norm(Z)
  % s.t.  S = Z

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

    % S subproblem. Linear least squares.
    % arg min_{S, Z} 1/2 sum_ti ||w_ti - R_t s_ti||^2 + rho/2 ||S - Z + U||_F^2
    V = Z - U;
    for t = 1:F
      R_t = rotations(:, :, t);
      W_t = projections(:, :, t);
      V_t = V(:, :, t);
      S_t = (R_t' * R_t + rho * eye(3)) \ (R_t' * W_t + rho * V_t);
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

    fprintf('%12d %12g %12g %12g\n', num_iter, rho, norm_r, norm_s);

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
