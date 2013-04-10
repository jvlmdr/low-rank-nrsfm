% Solves
%   min_x  ||P x - q||  s.t.  F x = g, X >= 0
% where x = vec(X) and X is symmetric (and square).

function X = semidefinite_constrained_least_squares(P, q, F, g)
  m = size(P, 1);
  n = sqrt(size(P, 2));
  k = size(F, 1);
  p = size(F, 2);

  % min_{x, e, r}  r  s.t. ||e|| <= r, P x - q = e, F x = g, X >= 0

  % min_{x, e, r}  r  s.t. [0, -I, P; 0, 0, F] [r; e; x] = [q; g], r >= ||e||,
  %                        X >= 0

  % Let z = [r; e; x]
  c_x = zeros(n * n, 1);
  c_e = zeros(m, 1);
  c_r = 1;
  c = [c_r; c_e; c_x];

  A_x = [P; F];
  A_e = [-eye(m); zeros(k, m)];
  A_r = zeros(m + k, 1);
  A = [A_r, A_e, A_x];
  b = [q; g];

  K = struct('q', [1 + m], 's', n);

  z = sedumi(A, b, c, K);

  x = z(m+2:end);
  X = mat(x, n);
end
