function S = find_structure_corrective_matrix(A_hat, B_hat, G1)

F = size(A_hat, 1) / 2;
K = size(B_hat, 1) / 3;
N = size(B_hat, 2);
n = 9 * K + F;
m = 6 * F;

A1 = A_hat * G1;

G = zeros(3 * K, 3, K);
C = zeros(F, K);
G(:, :, 1) = G1;
C(:, 1) = 1;

% Given G1, we need to find Gk for k = 2..K.
for k = 2:K
  % [ b_1 A_11 ]
  % [   ...    ] = A_hat G
  % [ b_F A_1F ]

  % Q c = L g

  % vec(kron(c, X)) = T c
  T = kron(ones(3, 1), kron(eye(F), ones(2, 1)));
  Q = diag(A1(:)) * T;

  L = kron(eye(3), A_hat);

  G_prev = reshape(G(:, :, 1:k-1), [3 * K, 3 * (k - 1)]);
  [U, D, V] = svd(G_prev, 'econ');
  E = kron(eye(3), U(:, 1:3*(k-1))');

  M = [                            Q, -L; ...
       zeros(size(E, 1), size(Q, 2)),  E];
  [U, D, V] = svd(M, 'econ');
  d = diag(D);
  dim = sum(d / d(1) < 1e-6);
  fprintf('dim(null(M)) = %d\n', dim);
  v = V(:, n);

  c = v(1:F);
  g = v(F + (1:9*K));
  Gk = reshape(g, [3 * K, 3]);

  C(:, k) = c;
  G(:, :, k) = Gk;
end

G = reshape(G, [3 * K, 3 * K]);
B = pinv(G) * B_hat;

S = kron(C, eye(3)) * B;

end
