% Solves
% arg min_{M, B} ||W - MB||_F  s.t.  M \in {motion matrices}
%
% [structure, rotations, basis, coeff] =
%   nrsfm_balm_metric_projections(projections, rotations, basis, max_iter)
%
% Parameters:
% projections -- 2 x P x F
% rotations -- 2 x 3 x F
% coeff -- K x F
% max_iter -- Number of iterations of alternation to do
%
% Returns:
% structure -- 3 x P x F
% rotations -- 2 x 3 x F
% basis -- 3 x K x P
% coeff -- K x F

function [structure, rotations, basis, coeff] = ...
    nrsfm_balm_metric_projections(projections, rotations, coeff, rho, ...
      max_iter, mu, tau_incr, tau_decr)
  % Introduce the auxiliary variable N.
  % arg min_{M, N, B} ||W - MB||_F s.t. N \in {motion matrices}, M - N = 0

  P = size(projections, 2);
  F = size(projections, 3);
  K = size(coeff, 1);

  % [K, F] -> [F, K]
  C = coeff';
  % [2, 3, F] -> [2F, 3F]
  R = block_diagonal_cameras(rotations);
  % Compose motion matrix.
  M = full(R * kron(C, eye(3)));
  % [2F, 3K] -> [2, F, 3, K]
  M = reshape(M, [2, F, 3, K]);

  W = projections_to_matrix(projections);

  converged = false;
  num_iter = 0;

  N = M;
  U = M - N;

  while ~converged && num_iter < max_iter
    prev_N = N;

    % (M, B) subproblem
    % arg min_{M, B} 1/2 ||W - MB||_F^2 + 1/2 rho ||M - N + U||_F^2
    V = N - U;
    V = reshape(V, [2 * F, 3 * K]);
    [M, B] = svd_proximity_operator(W, V, rho);
    M = reshape(M, [2, F, 3, K]);

    % N subproblem. "NRSFM Manifold Projector"
    % arg min_{N} ||M - N + U||_F^2  s.t.  N \in {motion matrices}
    V = M + U;
    % [2, F, 3, K] -> [2, 3, K, F]
    V = permute(V, [1, 3, 4, 2]);
    for t = 1:F
      V_t = V(:, :, :, t);
      [N_t, c_t, R_t] = project_motion_manifold(V_t);
      N(:, t, :, :) = N_t;
      % For output.
      coeff(:, t) = c_t;
      rotations(:, :, t) = R_t;
    end

    % Update multipliers and residuals.
    R = M - N;
    U = U + R;
    S = rho * (N - prev_N);

    norm_r = norm(R(:));
    norm_s = norm(S(:));

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

  basis = basis_from_matrix(B);
  %basis = find_basis(projections, rotations, coeff);

  structure = compose_structure(basis, coeff);
end
