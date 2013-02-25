% Returns a 2F x 3K(3K+1)/2 matrix of contrains.

function [A, c] = rotation_constraints(M_hat)
  % M_hat is 2F x 3K.
  F = size(M_hat, 1) / 2;
  K = size(M_hat, 2) / 3;

  % M = M_hat * G
  % M_t = M_hat_t G,  for all t

  % M_t M_t' = c_t^2 I
  % M_hat_t G G' M_hat_t' = c_t^2 I

  % Let Q = G G', Q = Q', Q is positive semidefinite.
  % M_hat_t Q M_hat_t' = c_t^2 I
  % kron(M_hat_t, M_hat_t) vec(Q) = vec(c_t^2 I)

  % [kron(M_hat_t, M_hat_t) vec(Q)]_1 - [kron(M_hat_t, M_hat_t) vec(Q)]_4 = 0
  % [kron(M_hat_t, M_hat_t) vec(Q)]_2 = 0
  % (Since the equations are symmetric, third constraint is redundant.)

  % Come back and enforce symmetry at the end.
  n = (3 * K) * (3 * K + 1) / 2;
  A = zeros(2, F, n);
  c = zeros(n, 1);

  % Make row pairs of M_hat easily accessible.
  M_hat = reshape(M_hat, [2, F, 3 * K]);
  % [2, F, 3K] -> [2, 3K, F]
  M_hat = permute(M_hat, [1, 3, 2]);

  H = construct_symmetric(3 * K);

  for t = 1:F
    A_t = kron(M_hat(:, :, t), M_hat(:, :, t));
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
