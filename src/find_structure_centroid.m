% Finds S that minimizes
%   || reshape(S) ||_{*}
% subject to
%   W = R (S + kron(1 1', mu))

function [X, mu] = find_structure_centroid(projections, use_3P, settings)

  N = length(projections.tracks);
  F = projections.num_frames;

  % Merge equations for all points.
  fprintf('Building linear system...\n');
  [A, b] = projection_equations_to_single_equation(projections);

  % Incorporate centroid.
  fprintf('Incorporating centroid...\n');
  C = [speye(3 * F * N), kron(ones(F * N, 1), speye(3))];
  A = A * C;
  n = 3 * F * N + 3;

  fprintf('Solving ADMM...\n');

  converged = false;
  num_iter = 0;
  rho = settings.rho;

  X = zeros(3 * F, N);
  Z = zeros(3 * F, N);
  U = X - Z;

  while ~converged && num_iter < settings.max_iter
    prev_Z = Z;

    % Constrained least-squares.
    V = Z - U;
    P = spdiags([ones(3 * F * N, 1); zeros(3, 1)], 0, n, n);
    q = [-V(:); zeros(3, 1)];
    m = size(A, 1);
    M = [P, A'; A, sparse(m, m)];
    x = M \ [-q; b];
    X = reshape(x(1:3 * F * N), [3 * F, N]);
    mu = x(3 * F * N + (1:3));

    % Singular value soft thresholding.
    V = X + U;
    if use_3P
      V = reshape(permute(reshape(V, [3, F, N]), [1, 3, 2]), [3 * N, F]);
    end
    Z = singular_value_soft_threshold(V, 1 / rho);
    if use_3P
      Z = reshape(permute(reshape(Z, [3, N, F]), [1, 3, 2]), [3 * F, N]);
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

  X = X + kron(ones(F, N), mu);
end
