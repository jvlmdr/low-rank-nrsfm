% Returns a m x 3K(3K+1)/2 system of orthogonality and basis constraints
% with m = F(K-1) - K(K-1)/2 = (2F-K)(K-1)/2 the number of constraints.

function [A, b] = rotation_and_basis_constraints(M_hat, k)
  % M_hat is 2F x 3K.
  F = size(M_hat, 1) / 2;
  K = size(M_hat, 2) / 3;

  % Number of variables.
  n = 3 * K * (3 * K + 1) / 2;

  % Make row pairs of M_hat easily accessible.
  % [2F, 3K] -> [2, F, 3K] -> [2, 3K, F]
  M_hat = reshape(M_hat, [2, F, 3 * K]);
  M_hat = permute(M_hat, [1, 3, 2]);

  H = construct_symmetric(3 * K);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Pairs (t, t) such that t <= K.

  A_inside = zeros(3, n, K);
  b_inside = zeros(3, K);

  for t = 1:K
    % [M_hat_t Q M_hat_t'] = c_tk I, with c_tk = delta[t - k].
    c_tk = double(t == k);
    % LHS is symmetric.

    % [kron(M_hat_t, M_hat_t) vec(Q)]_1 = c_tk
    % [kron(M_hat_t, M_hat_t) vec(Q)]_2 = 0
    % [kron(M_hat_t, M_hat_t) vec(Q)]_4 = c_tk
    A_t = kron(M_hat(:, :, t), M_hat(:, :, t));
    % Convert to symmetric parametrization.
    A_t = A_t * H;

    % Append to system.
    A_inside(:, :, t) = [A_t(1, :); A_t(2, :); A_t(4, :)];
    b_inside(:, t) = [c_tk; 0; c_tk];
  end

  % [3, n, K] -> [3, K, n] -> [3K, n]
  A_inside = permute(A_inside, [1, 3, 2]);
  A_inside = reshape(A_inside, [3 * K, n]);
  b_inside = b_inside(:);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Pairs (t, t) such that t > K.

  A_outside = zeros(2, n, F - K);
  b_outside = zeros(2, F - K);

  for t = (K + 1):F
    % [M_hat_t Q M_hat_t'] = c_tk I, with c_tk unknown.
    % LHS is symmetric.

    % [kron(M_hat_t, M_hat_t) vec(Q)]_1 - [kron(M_hat_t, M_hat_t) vec(Q)]_4 = 0
    % [kron(M_hat_t, M_hat_t) vec(Q)]_2 = 0
    A_t = kron(M_hat(:, :, t), M_hat(:, :, t));
    % Convert to symmetric parametrization.
    A_t = A_t * H;

    % Append to system.
    A_outside(:, :, t - K) = [A_t(1, :) - A_t(4, :); A_t(2, :)];
    b_outside(:, t - K) = [0; 0];
  end

  % [2, n, F-K] -> [2, F-K, n] -> [2(F-K), n]
  A_outside = permute(A_outside, [1, 3, 2]);
  A_outside = reshape(A_outside, [2 * (F - K), n]);
  b_outside = b_outside(:);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Unique unordered pairs (t, u), t != u, such that at least one of t, u
  % satisfies v <= K, v != k.

  m = (K - 1) * K / 2 + (K - 1) * (F - K);
  A_cross = zeros(4, n, m);
  b_cross = zeros(4, m);
  i = 1;

  for t = [1:(k - 1), (k + 1):K]
    if t < k
      range = [(t + 1):k, (K + 1):F];
    else
      range = [1:(t - 1), (K + 1):F];
    end

    for u = range
      % [M_hat_t Q M_hat_u'] = 0
      % LHS is not symmetric.

      % [kron(M_hat_u, M_hat_t) vec(Q)] = 0
      A_tu = kron(M_hat(:, :, u), M_hat(:, :, t));
      % Convert to symmetric parametrization.
      A_tu = A_tu * H;

      % Append to system.
      A_cross(:, :, i) = A_tu;
      b_cross(:, i) = zeros(4, 1);
      i = i + 1;
    end
  end

  % [4, n, m] -> [4, m, n] -> [4m, n]
  A_cross = permute(A_cross, [1, 3, 2]);
  A_cross = reshape(A_cross, [4 * m, n]);
  b_cross = b_cross(:);

  A = [A_inside; A_outside; A_cross];
  b = [b_inside; b_outside; b_cross];
end
