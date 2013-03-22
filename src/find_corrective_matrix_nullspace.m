% Recovers corrective matrix and basis coefficients given M_hat and R.
%
% [G, C] = find_corrective_matrix_nullspace(M_hat, R)
%
% Parameters:
% M_hat -- 2F x 3K matrix from SVD
% R -- 2F x 3F

function [G, C] = find_corrective_matrix_nullspace(M_hat, R)
  F = size(M_hat, 1) / 2;
  K = size(M_hat, 2) / 3;

  % Find G and C that minimize || M_hat G - R kron(C, I) ||

  % (R kron(C, I))'
  %   = kron(C', I) R'
  %   = kron(C', I) [r1 ... rn]
  %   = [kron(C', I) r1 ... kron(C', I) rn]
  %   = [vec(R1 C) ... vec(Rn C)]
  %   = [kron(I, R1) c ... kron(I, Rn) c]
  % vec[(R kron(C, I))']
  %   = [kron(I, R1) c; ...; kron(I, Rn) c]
  %   = [kron(I, R1); ...; kron(I, Rn)] c

  % Extract rows of R.
  r = mat2cell(R, ones(2 * F, 1), 3 * F);
  % Convert from 3F vectors to 3xF matrices.
  r = cellfun(@(x) { reshape(x, [3, F]) }, r);
  % Build full system.
  A_c = cell2mat(r);

  % Re-order rows for transpose.
  order = 1:(6 * F);
  order = reshape(reshape(order, [3, 2 * F])', [6 * F, 1]);
  A_c = A_c(order, :);

  % Build matrix for G.
  A_q = kron(speye(3), M_hat);

  % Assemble full system.
  A = [A_q, -A_c];

  [U, S, V] = svd(full(A));

  % Take last K singular vectors.
  V = V(:, (end - K + 1):end);
  G = V(1:(9 * K), :);
  C = V((9 * K + 1):end, :);

  G = reshape(G, [3 * K, 3 * K]);
end
