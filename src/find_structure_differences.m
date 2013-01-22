% Finds S that minimizes
%   || reshape(D3 * S) ||_{*}
% subject to
%   W = R * S
% where D3 = kron(D, eye(3)) and D is the (F-1) x F first-difference matrix

function X = find_structure_differences(projections, use_3P, settings)

  N = length(projections.tracks);
  F = projections.num_frames;

  fprintf('Building linear systems...\n');
  systems = projection_equations_to_trajectory_equations(projections);

  D = kron(first_difference_matrix(F), speye(3));
  P = D' * D;

  converged = false;
  num_iter = 0;
  rho = settings.rho;

  X = zeros(3 * F, N);
  Z = zeros(3 * (F - 1), N);
  U = D * X - Z;

  while ~converged && num_iter < settings.max_iter
    prev_Z = Z;

    % Constrained least-squares, independent per point.
    V = Z - U;
    for i = 1:N
      q = -D' * V(:, i);
      A = systems(i).A;
      b = systems(i).b;
      m = size(A, 1);
      M = [P, A'; A, sparse(m, m)];
      x = M \ [-q; b];
      X(:, i) = x(1:3 * F);
    end

    % Singular value soft thresholding.
    V = D * X + U;
    if use_3P
      V = reshape(permute(reshape(V, [3, F - 1, N]), [1, 3, 2]), ...
          [3 * N, F - 1]);
    end
    Z = singular_value_soft_threshold(V, 1 / rho);
    if use_3P
      Z = reshape(permute(reshape(Z, [3, N, F - 1]), [1, 3, 2]), ...
          [3 * (F - 1), N]);
    end

    % Update multipliers.
    R = D * X - Z;
    U = U + R;

    S = rho * D' * (Z - prev_Z);

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
end
