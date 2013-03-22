% Finds the corrective matrix and basis coefficients as the solution to a
% homogeneous linear system given cameras and an estimate of the basis.
%
% [G, C] = find_corrective_matrix(M_hat, R, B)
%
% Parameters:
% M_hat -- 2F x 3K matrix from SVD
% R -- 2F x 3F
% B -- 3K x P

function [G, C] = find_corrective_matrix(M_hat, R, B)
  F = size(M_hat, 1) / 2;
  K = size(M_hat, 2) / 3;
  P = size(B, 2);

  % Find G and C that minimize || M_hat G B - R kron(C, I) B ||_F

  % (R kron(C, I))'
  %   = kron(C', I) R'
  %   = kron(C', I) [r1 ... rm]
  %   = [kron(C', I) r1 ... kron(C', I) rm]
  %   = [vec(R1 C) ... vec(Rm C)]
  %   = [kron(I, R1) c ... kron(I, Rm) c]

  % vec((R kron(C, I))')
  %   = [kron(I, R1); ...; kron(I, Rm)] c

  % vec((R kron(C, I) B)')
  %   = vec(B' (R kron(C, I))')
  %   = kron(I, B') vec((R kron(C, I))')
  %   = kron(I, B') [kron(I, R1); ...; kron(I, Rm)] c

  % Extract rows of R.
  r = mat2cell(R, ones(2 * F, 1), 3 * F);
  % Convert from 3F vectors to 3xF matrices.
  r = cellfun(@(x) { kron(speye(K), reshape(x, [3, F])) }, r);
  % Build full system.
  A_c1 = cell2mat(r);
  A_c2 = kron(speye(2 * F), B');
  A_c = A_c2 * A_c1;

  % Build matrix for G.
  % vec((M_hat G B)') = vec(B' G' M_hat') = kron(M_hat, B') vec(G')
  A_g = kron(M_hat, B');

  % Assemble full system.
  A = [A_g, -A_c];

  g_begin = 1;
  g_end = g_begin + 9 * K ^ 2;
  c_begin = g_end;
  c_end = c_begin + F * K;

  % Get last singular vector.
  n = size(A, 2);
  opts = struct('issym', true, 'isreal', true, 'disp', 0);
  [V, D] = eigs(@(x) A' * (A * x), n, 1, 'SA', opts);
  %AA = A' * A;
  %[V, D] = eigs(AA, 1, 'SA');

  v = V(:, end);
  g = v(g_begin:g_end - 1);
  c = v(c_begin:c_end - 1);
  G = reshape(g, [3 * K, 3 * K])';
  C = reshape(c, [F, K]);

  fprintf('norm(A v): %g\n', norm(A * v));
  fprintf('cond(G): %g\n', cond(G));
end
