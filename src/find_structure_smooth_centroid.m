% Finds S that minimizes
%   || reshape(S) ||_{*} + 1/2 || D S ||_{F}^{2}
% subject to
%   W = R S

function [X, mu] = find_structure_smooth_centroid(projections, use_3P, D, ...
    lambda, settings)

  N = length(projections);
  F = size(projections(1).Q, 2) / 3;

  % Merge equations for all points and centroid.
  n = 3 * F * N + 3;
  A = sparse(0, n);
  b = zeros(0, 1);
  for i = 1:N
    Q = projections(i).Q;
    % Add centroid mu in every frame.
    A_mu = kron(ones(F, 1), speye(3));
    A_mu = Q * A_mu;
    m = size(Q, 1);
    A_x = sparse(m, 3 * F * N);
    A_x(:, (3 * F * (i - 1) + 1):(3 * F * i)) = Q;
    A = [A; A_x, A_mu];
    b = [b; projections(i).q];
  end

  D = kron(speye(N), kron(D, speye(3)));
  DD = D' * D;

  converged = false;
  num_iter = 0;
  rho = settings.rho;

  X = zeros(3 * F, N);
  Z = zeros(3 * F, N);
  U = X - Z;

  while ~converged && num_iter < settings.max_iter
    num_iter = num_iter + 1;
    prev_Z = Z;

    % Constrained least-squares.
    V = Z - U;
    P_X = rho * speye(3 * F * N) + lambda * DD;
    P = [P_X, sparse(3 * F * N, 3); sparse(3, 3 * F * N), sparse(3, 3)];
    q = [-rho * V(:); zeros(3, 1)];
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

    if norm_r ~= 0 && norm_s ~= 0
      if norm_r > settings.mu * norm_s
        rho = rho * settings.tau_incr;
      elseif norm_s > settings.mu * norm_r
        rho = rho / settings.tau_decr;
      end
    end
  end

  X = X + kron(ones(F, N), mu);
end
