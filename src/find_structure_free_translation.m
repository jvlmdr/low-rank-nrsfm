% Finds S that minimizes
%   || reshape(S) ||_{*}
% subject to
%   W = R (S + mu * 1')
%
% where mu is a 3Fx1 centroid trajectory.
%
% Parameters:
% projections
% - num_frames
% - tracks
%   - frames
%   - equations
%     - A, b -- Linear system

function [X, mu] = find_structure_free_translation(projections, use_3P, ...
    settings)

  N = length(projections.tracks);
  F = projections.num_frames;

  systems = projection_equations_to_shape_equations(projections);

  % Incorporate centroid into each system.
  C = [speye(3 * N), kron(ones(N, 1), speye(3))];
  for t = 1:F
    systems(t).A = systems(t).A * C;
  end

  converged = false;
  num_iter = 0;
  rho = settings.rho;

  X = zeros(3 * N, F);
  Z = randn(3 * N, F);
  U = X - Z;
  mu = zeros(3, F);

  while ~converged && num_iter < settings.max_iter
    prev_Z = Z;

    % Constrained least-squares subproblem.
    V = Z - U;
    for t = 1:F
      n = 3 * N + 3;
      P = spdiags([ones(3 * N, 1); zeros(3, 1)], 0, n, n);
      q = [-V(:, t); zeros(3, 1)];
      A = systems(t).A;
      b = systems(t).b;

      m = size(A, 1);
      M = [P, A'; A, sparse(m, m)];
      x = M \ [-q; b];

      X(:, t) = x(1:3 * N);
      mu(:, t) = x(3 * N + (1:3));
    end

    % Singular value soft thresholding.
    V = X + U;
    if ~use_3P
      V = reshape(permute(reshape(V, [3, N, F]), [1, 3, 2]), [3 * F, N]);
    end
    Z = singular_value_soft_threshold(V, 1 / rho);
    if ~use_3P
      Z = reshape(permute(reshape(Z, [3, F, N]), [1, 3, 2]), [3 * N, F]);
    end

    % Update multipliers.
    R = X - Z;
    U = U + R;

    S = rho * (Z - prev_Z);

    norm_r = norm(R(:));
    norm_s = norm(S(:));

    fprintf('%12d %12g %12g %12g\n', num_iter, rho, norm_r, norm_s);

    if num_iter > 0 && num_iter < settings.max_iter / 2
      if norm_r ~= 0 && norm_s ~= 0
        if norm_r > settings.mu * norm_s
          rho = rho * settings.tau_incr;
        elseif norm_s > settings.mu * norm_r
          rho = rho / settings.tau_decr;
        end
      end
    end

    num_iter = num_iter + 1;
  end

  X = reshape(permute(reshape(X, [3, N, F]), [1, 3, 2]), [3 * F, N]);
  mu = mu(:);
  X = X + mu * ones(1, N);
end
