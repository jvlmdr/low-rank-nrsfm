function G = find_rotations_using_first_triple(P_hat, lambda)
  % Find G = Q Q', where Q is the 3K x 3 matrix which best corrects P_hat.

  % ADMM settings.
  settings = struct(...
      'rho', 1, ...
      'mu', 10, ...
      'tau_incr', 2, ...
      'tau_decr', 2, ...
      'max_iter', 80, ...
      'epsilon_abs', 1e-3, ...
      'epsilon_rel', 1e-3, ...
      'min_rho_iter', 4);

  F = size(P_hat, 1) / 2;
  K = size(P_hat, 2) / 3;

  n = (3 * K) * (3 * K + 1) / 2;

  % Get constraints that make A orthogonal.
  A = rotation_constraints(P_hat);
  m = size(A, 2);
  % c' x = d.
  c = zeros(n, 1);
  c(1) = 1;
  C = c';
  d = 1;

  AA = A' * A;

  % Maps a vector of size n to a vectorized matrix of size 3K x 3K.
  H = construct_symmetric(3 * K);

  % trace(X) = g' vec(X)
  g = (ones(3 * K, 1)' * diag_vec(3 * K) * H)';

  z = zeros(n, 1);
  u = zeros(n, 1);

  converged = false;
  num_iter = 0;
  rho = settings.rho;
  num_rho_iter = 0;

  while ~converged && num_iter < settings.max_iter
    z_prev = z;

    % Solve x subproblem:
    v = z - u;
    % System of equations with trivial solution.
    P = lambda * AA + rho * eye(n);
    q = g - rho * v;
    x = P \ -q;

    % Solve z subproblem:
    v = x + u;
    V = reshape(H * v, [3 * K, 3 * K]);
    Z = project_psd(V);
    z = H \ Z(:);

    r = x - z;
    u = u + r;
    s = rho * (z - z_prev);

    norm_r = norm(r);
    norm_s = norm(s);

    fprintf('%12d %12g %12g %12g\n', num_iter, rho, norm_r, norm_s);

    num_iter = num_iter + 1;
    num_rho_iter = num_rho_iter + 1;

    if ~(num_rho_iter < settings.min_rho_iter) && ...
        num_iter < settings.max_iter / 2
      if norm_r ~= 0 && norm_s ~= 0
        if norm_r > settings.mu * norm_s
          rho = rho * settings.tau_incr;
          % Keep y = rho * u constant.
          u = u / settings.tau_incr;
          num_rho_iter = 0;
        elseif norm_s > settings.mu * norm_r
          rho = rho / settings.tau_decr;
          % Keep y = rho * u constant.
          u = u * settings.tau_decr;
          num_rho_iter = 0;
        end
      end
    end
  end

  G = reshape(H * z, [3 * K, 3 * K]);
end
