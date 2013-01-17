% Finds the min nuclear norm F x 3P structure
% Projections is enforced by W = R S

function S = find_structure_shape_basis_admm(W, R)

rho_init = 1;
mu = 10;
tau_incr = 2;
tau_decr = 2;
max_iter = 200;

F = size(W, 1) / 2;
N = size(W, 2);
n = 3 * F * N;
m = size(R, 1) * N;

order = 1:n;
order = reshape(order, [3, F, N]);
order = permute(order, [2, 1, 3]);
order = order(:);

B = speye(n, n);
B = B(order, :);

converged = false;
num_iter = 0;

A = kron(speye(N), R);
b = W(:);

rho = rho_init;
z = zeros(n, 1);
u = zeros(n, 1);

while ~converged && num_iter < max_iter
  z_prev = z;
  num_iter = num_iter + 1;

  % Solve x subproblem.
  v = B * z - u;
  V = reshape(v, [F, 3 * N]);
  X = singular_value_soft_threshold(V, 1 / rho);
  stem(svd(X));
  drawnow;
  x = X(:);

  % Solve z subproblem.
  v = B' * (x + u);
  P = rho * speye(n);
  q = -rho * v;
  zero = sparse(m, m);
  s = [P, A'; A, zero] \ [-q; b];
  z = s(1:n);

  r = x - B * z;
  u = u + r;
  s = rho * -B * (z - z_prev);

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

S = reshape(B' * x, [3 * F, N]);
%S = reshape(z, [3 * F, N]);

end
