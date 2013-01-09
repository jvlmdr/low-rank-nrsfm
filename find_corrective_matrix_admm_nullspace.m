function G = find_corrective_matrix_admm_nullspace(PI_Hat)

epsilon = 1e-6;
rho_init = 1;
mu = 10;
tau_incr = 2;
tau_decr = 2;

% PI_Hat is 2F x 3K.
F = size(PI_Hat, 1) / 2;
K = size(PI_Hat, 2) / 3;

% Solving for X which is n x n.
n = 3 * K;
% Unique entries in symmetric X.
m = n * (n + 1) / 2;

[A_rot, ~, A_scale] = rotation_constraints(PI_Hat);
A_rot = A_rot * construct_symmetric(n);
b_rot = zeros(size(A_rot, 1), 1);
A_scale = A_scale * construct_symmetric(n);
b_scale = 2 * m;

% Additional constraint that sum(i' i + j' j) = 1.
A = [A_rot; A_scale];
b = [b_rot; b_scale];

c = -(A \ b);

[U, D, V] = svd(A);
d = diag(D);
phi = V(:, d / d(1) < epsilon);
k = size(phi, 2);
fprintf('"Nullspace" has %d dimensions\n', k);

% trace(X) = c' vec(X)
g = (ones(n, 1)' * diag_vec(n) * construct_symmetric(n))';

z = zeros(m, 1);
u = zeros(m, 1);
rho = rho_init;

converged = false;
num_iter = 0;

while ~converged && num_iter < 400
  num_iter = num_iter + 1;
  z_prev = z;

  % Solve x subproblem:
  v = z + c - u;

  P = rho * eye(k);
  q = phi' * (g - rho * v);

  alpha = P \ -q;

  % Solve z subproblem:
  v = phi * alpha - c + u;
  V = reshape(construct_symmetric(n) * v, [n, n]);
  Z = project_psd(V);
  plot(svd(Z));
  drawnow;
  z = construct_symmetric(n) \ Z(:);

  r = phi * alpha - z - c;
  u = u + r;
  s = rho * phi' * -eye(m) * (z - z_prev);

  norm_r = norm(r);
  norm_s = norm(s);

  fprintf('%12d %12g %12g %12g\n', num_iter, rho, norm_r, norm_s);

  if norm_r ~= 0 && norm_s ~= 0
    if norm_r > mu * norm_s
      rho = rho * tau_incr;
    elseif norm_s > mu * norm_r
      rho = rho / tau_decr;
    end
  end
end

x = phi * alpha;
X = reshape(construct_symmetric(n) * x, [n, n]);

[U_Hat D_Hat V_Hat] = svd(X);

G_Hat = U_Hat(:,1:3)*sqrt(D_Hat(1:3,1:3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nonlinear optimization
disp('  Nonlinear refinement of G...(wait).'); 

options = optimset('Display', 'Final', 'Diagnostics','off','Largescale', 'off', 'MaxFunEval',200000,'MaxIter',5000,'TolFun',1e-10,'TolX',1e-10);

[G, fval] = fminunc(@evalQ_regularization,G_Hat,options,PI_Hat);  %fminunc  fminsearch

end
