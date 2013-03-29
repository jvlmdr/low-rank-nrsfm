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

%%  cvx_solver sedumi;
%%  cvx_begin sdp quiet;
%%    variable X(6, 6) symmetric
%%    A = X(1:3, 1:3);
%%    B = X(1:3, 4:6);
%%    C = X(4:6, 4:6);
%%    w = [B(2, 3) - B(3, 2); B(3, 1) - B(1, 3); B(1, 2) - B(2, 1)];
%%
%%    minimize trace(E * X)
%%    subject to
%%      X >= 0
%%      trace(A) == 1;
%%      trace(C) == 1;
%%      trace(B) == 0;
%%      [eye(3) - A - C, w; w', 1] >= 0;
%%  cvx_end;

  % One 6x6 and one 4x4 symmetric positive semidefinite matrix.
  % X = [A B; B' C]
  % Y = [I - A - C, w; w' 1]
  % w = [B(2, 3) - B(3, 2); B(3, 1) - B(1, 3); B(1, 2) - B(2, 1)];
  cones = struct('s', [6, 4]);

  % Objective.
  c = [vec(E); zeros(4 * 4, 1)];

  % Linear equality constraints.
  % trace(A) == 1
  a1 = vec(diag([1, 1, 1, 0, 0, 0]));
  b1 = 1;
  % trace(C) == 1
  a2 = vec(diag([0, 0, 0, 1, 1, 1]));
  b2 = 1;
  % trace(B) == 0
  a3 = vec(diag([1, 1, 1], 3));
  b3 = 0;

  A_tr = [[a1'; a2'; a3'], zeros(3, 4 * 4)];
  b_tr = [b1; b2; b3];

  % Y = [I - A - C, w; w' 1]
  % Y + [A + C, -w] = [I, 0]
  %     [  -w',  0]   [0, 1]
  %
  % Y = [Q, q; q', p]
  % [Q + A + C, q - w] = [I, 0]
  % [ (q - w)',     p]   [0, 1]
  % Extract top-left 3x3 matrix of 4x4 matrix.
  Q_Y = zeros(3, 3, 4, 4);
  Q_Y(:, :, 1:3, 1:3) = reshape(eye(3 * 3), [3, 3, 3, 3]);
  % Extract top-left 3x3 matrix of 6x6 matrix.
  A_X = zeros(3, 3, 6, 6);
  A_X(:, :, 1:3, 1:3) = reshape(eye(3 * 3, 3 * 3), [3, 3, 3, 3]);
  % Extract bottom-right 3x3 matrix of 6x6 matrix.
  C_X = zeros(3, 3, 6, 6);
  C_X(:, :, 4:6, 4:6) = reshape(eye(3 * 3, 3 * 3), [3, 3, 3, 3]);
  % Vectorize.
  A_X = reshape(A_X, [3 * 3, 6 * 6]);
  C_X = reshape(C_X, [3 * 3, 6 * 6]);
  Q_Y = reshape(Q_Y, [3 * 3, 4 * 4]);

  % Construct w from vec(B).
  % w = [B(2, 3) - B(3, 2); B(3, 1) - B(1, 3); B(1, 2) - B(2, 1)];
  w_B = zeros(3, 3, 3);
  w_B(1, 2, 3) = 1;
  w_B(1, 3, 2) = -1;
  w_B(2, 3, 1) = 1;
  w_B(2, 1, 3) = -1;
  w_B(3, 1, 2) = 1;
  w_B(3, 2, 1) = -1;
  % Construct w from X.
  w_X = zeros(3, 6, 6);
  w_X(:, 1:3, 4:6) = w_B;
  % Extract top-right 3x1 vector of Y.
  q_Y = zeros(3, 4, 4);
  q_Y(:, 1:3, 4) = eye(3);
  % Vectorize.
  w_X = reshape(w_X, [3, 6 * 6]);
  q_Y = reshape(q_Y, [3, 4 * 4]);

  A = [          A_tr;
       A_X + C_X, Q_Y;
            -w_X, q_Y];
  b = [       b_tr;
       vec(eye(3));
       zeros(3, 1)];

  pars = struct('fid', 0);

  x = sedumi(A, b, c, cones, pars);

  X = x(1:(6 * 6));
  X = reshape(X, [6, 6]);

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
