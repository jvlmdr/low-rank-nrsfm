function [A, c] = rotation_constraints(P_hat)
  % P_hat is 2F x 3K.
  F = size(P_hat, 1) / 2;
  K = size(P_hat, 2) / 3;

  % P = P_hat * Q
  % P_t = P_hat_t Q,  for all t

  % P_t P_t' = c_t^2 I
  % P_hat_t Q Q' P_hat_t' = c_t^2 I

  % Let G = Q Q', G = G', G is positive semidefinite.
  % P_hat_t G P_hat_t' = c_t^2 I
  % kron(P_hat_t, P_hat_t) vec(G) = vec(c_t^2 I)

  % [kron(P_hat_t, P_hat_t) vec(G)]_1 - [kron(P_hat_t, P_hat_t) vec(G)]_4 = 0
  % [kron(P_hat_t, P_hat_t) vec(G)]_2 = 0
  % (Since the equations are symmetric, third constraint is redundant.)

  % Come back and enforce symmetry at the end.
  n = (3 * K) * (3 * K + 1) / 2;
  A = zeros(2, F, n);
  c = zeros(n, 1);

  P_hat = reshape(P_hat, [2, F, 3 * K]);
  P_hat = permute(P_hat, [1, 3, 2]);

  H = construct_symmetric(3 * K);

  for t = 1:F
    A_t = kron(P_hat(:, :, t), P_hat(:, :, t));
    % Convert to symmetric.
    A_t = A_t * H;
    % Append to system.
    A(:, t, :) = [A_t(1, :) - A_t(4, :); A_t(2, :)];

    % Magnitude constraint.
    if t == 1
      c = A_t(1, :)';
    end
  end
  A = reshape(A, [2 * F, n]);
end
