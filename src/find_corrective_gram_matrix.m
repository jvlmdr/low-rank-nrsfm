% Solves
%   min trace(X) + lambda || A x ||
%   s.t. X is positive semidefinite,
%        C x >= 1
% where X = symmetric_matrix(x).
%
% Note that X is p x p, A is m x n, with n = p (p + 1) / 2.

function Q = find_corrective_gram_matrix(A, C, lambda)
  m = size(A, 1);
  n = size(A, 2);
  % p (p + 1) = 2n
  % p < sqrt(2n) < p + 1
  p = floor(sqrt(2 * n));

  d = 1;

  % Maps a vector of size n to a vectorized matrix of size 3K x 3K.
  [H, subset] = construct_symmetric(p);

  % Solve using CVX.
  cvx_begin sdp quiet;
    cvx_solver sedumi;
    variable Q(p, p) symmetric;
    q = Q(subset)';

    minimize trace(Q) + lambda * norm(A * q);
    subject to
      Q >= 0;
      %C(1, :) * q == d;
      C * q >= d;
  cvx_end;

%  % ADMM settings.
%  settings = struct(...
%      'rho', 1, ...
%      'mu', 10, ...
%      'tau_incr', 2, ...
%      'tau_decr', 2, ...
%      'max_iter', 80, ...
%      'epsilon_abs', 1e-3, ...
%      'epsilon_rel', 1e-3, ...
%      'min_rho_iter', 4);
%
%  AA = A' * A;
%
%  % trace(X) = g' vec(X)
%  g = (ones(p, 1)' * diag_vec(p) * H)';
%
%  z = zeros(n, 1);
%  u = zeros(n, 1);
%
%  converged = false;
%  num_iter = 0;
%  rho = settings.rho;
%  num_rho_iter = 0;
%
%  while ~converged && num_iter < settings.max_iter
%    z_prev = z;
%
%    % Solve x subproblem:
%    v = z - u;
%    % TODO: Possible to avoid building Gram matrix?
%    P = lambda * AA + rho * eye(n);
%    q = g - rho * v;
%
%    s = [P, C'; C, 0] \ [-q; d];
%    x = s(1:n);
%
%    % Solve z subproblem:
%    v = x + u;
%    V = reshape(H * v, [p, p]);
%    Z = project_psd(V);
%    z = H \ Z(:);
%
%    r = x - z;
%    u = u + r;
%    s = rho * (z - z_prev);
%
%    norm_r = norm(r);
%    norm_s = norm(s);
%
%    fprintf('%12d %12g %12g %12g\n', num_iter, rho, norm_r, norm_s);
%
%    num_iter = num_iter + 1;
%    num_rho_iter = num_rho_iter + 1;
%
%    if ~(num_rho_iter < settings.min_rho_iter) && ...
%        num_iter < settings.max_iter / 2
%      if norm_r ~= 0 && norm_s ~= 0
%        if norm_r > settings.mu * norm_s
%          rho = rho * settings.tau_incr;
%          % Keep y = rho * u constant.
%          u = u / settings.tau_incr;
%          num_rho_iter = 0;
%        elseif norm_s > settings.mu * norm_r
%          rho = rho / settings.tau_decr;
%          % Keep y = rho * u constant.
%          u = u * settings.tau_decr;
%          num_rho_iter = 0;
%        end
%      end
%    end
%  end
%
%  Q = reshape(H * z, [p, p]);
end
