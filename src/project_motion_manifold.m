% Minimizes ||M - kron(c', R)||_F  s.t.  R R' = I.
%
% Parameters:
% M -- 2 x 3 x K
%
% Returns:
% M -- 2 x 3 x K
% c -- K x 1
% R -- 2 x 3

function [M, c, R] = project_motion_manifold(M)
  K = size(M, 3);

  E = zeros(6, 6);
  for i = 1:K
    M_i = M(:, :, i)';
    E = E - M_i(:) * M_i(:)';
  end

  cvx_begin sdp quiet;
    variable X(6, 6) symmetric
    A = X(1:3, 1:3);
    B = X(1:3, 4:6);
    C = X(4:6, 4:6);
    w = [B(2, 3) - B(3, 2); B(3, 1) - B(1, 3); B(1, 2) - B(2, 1)];

    minimize trace(E * X)
    subject to
      X >= 0
      trace(A) == 1;
      trace(C) == 1;
      trace(B) == 0;
      [eye(3) - A - C, w; w', 1] >= 0;
  cvx_end;

  [U, S, V] = svd(X, 'econ');
  q = sqrt(2) * V(:, 1);
  R = reshape(q, [3, 2])';

  % Find nearest rotation matrix.
  R = nearest_scaled_rotation_matrix(R);

  c = zeros(K, 1);
  for i = 1:K
    c(i) = 1/2 * trace(R * M(:, :, i)');
  end

  M_mat = kron(c', R);

  % [2, 3K] -> [2, 3, K]
  M = reshape(M_mat, [2, 3, K]);
end
