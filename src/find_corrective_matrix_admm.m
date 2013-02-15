function G = find_corrective_matrix_admm(pi_hat)
  lambda = 1e9;

  % ADMM settings.
  settings = struct(...
      'rho', 1, ...
      'mu', 10, ...
      'tau_incr', 2, ...
      'tau_decr', 2, ...
      'max_iter', 200, ...
      'epsilon_abs', 1e-3, ...
      'epsilon_rel', 1e-3, ...
      'min_rho_iter', 4);

  % pi_hat is 2F x 3K.
  F = size(pi_hat, 1) / 2;
  K = size(pi_hat, 2) / 3;

  % Solving for X which is n x n.
  n = 3 * K;
  % Unique entries in symmetric X.
  m = n * (n + 1) / 2;

  [C, ~, A] = rotation_constraints_dai(pi_hat);
  C = C * construct_symmetric(n);

  % Additional constraint that sum(i_t' i_t + j_t' j_t) = 2 F.
  A = A * construct_symmetric(n);
  b = 2 * F;

  % trace(X) = c' vec(X)
  c = (ones(n, 1)' * diag_vec(n) * construct_symmetric(n))';

  z = zeros(m, 1);
  u = zeros(m, 1);

  converged = false;
  num_iter = 0;
  rho = settings.rho;
  num_rho_iter = 0;

  while ~converged && num_iter < settings.max_iter
    z_prev = z;

    % Solve x subproblem:
    v = z - u;
    % TODO: Possible to avoid building Gram matrix?
    P = lambda * C' * C + rho * eye(m);
    q = c - rho * v;

    s = [P, A'; A, 0] \ [-q; b];
    x = s(1:m);

    % Solve z subproblem:
    v = x + u;
    V = reshape(construct_symmetric(n) * v, [n, n]);
    Z = project_psd(V);
    stem(svd(Z));
    drawnow;
    z = construct_symmetric(n) \ Z(:);

    r = x - z;
    u = u + r;
    s = rho * eye(m)' * -eye(m) * (z - z_prev);

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

  X = reshape(construct_symmetric(n) * x, [n, n]);

  [U_Hat D_Hat V_Hat] = svd(X);

  G_Hat = U_Hat(:,1:3)*sqrt(D_Hat(1:3,1:3));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Nonlinear optimization
  disp('  Nonlinear refinement of G...(wait).'); 

  options = optimset('Display', 'Final', 'Diagnostics','off','Largescale', 'off', 'MaxFunEval',200000,'MaxIter',5000,'TolFun',1e-10,'TolX',1e-10);

  [G, fval] = fminunc(@evalQ_regularization,G_Hat,options,pi_hat);  %fminunc  fminsearch
end
