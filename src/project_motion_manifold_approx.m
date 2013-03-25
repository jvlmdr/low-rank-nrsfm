% Approximately minimizes ||M - kron(c', R)||_F  s.t.  R R' = I.
%
% Parameters:
% M -- 2 x 3 x K
%
% Returns:
% M -- 2 x 3 x K
% c -- K x 1
% R -- 2 x 3

function [M, c, R] = project_motion_manifold_approx(M)
  K = size(M, 3);
  M_mat = reshape(M, [2, 3 * K]);

  % Take sum of Gram matrices.
  MM = zeros(3, 3);
  for i = 1:K
    MM = MM + M(:, :, i)' * M(:, :, i);
  end

  [V, D] = eig(MM);
  % Take two largest magnitude eigenvectors.
  d = diag(D);
  [~, ind] = sort(abs(d), 'descend');
  ind = ind(1:2);
  d = d(ind);
  V = V(:, ind)';

  % Now R = Q V, where V is 2 x 3 and Q is 2 x 2.
  R = zeros(2, 3, 2);
  c = zeros(K, 2);
  y = zeros(2, 1);

  for j = 1:2
    sgn = (-1) ^ j;
    F1 = [1,    0;
          0,  sgn];
    F2 = [0, -sgn;
          1,    0];

    AA = zeros(2, 2);

    for i = 1:K
      T = M(:, :, i) * V';
      A = T(:)' * [F1(:), F2(:)];

      AA = AA + A' * A;
    end

    [X, D] = eig(AA);
    d = abs(diag(D));
    x = X(:, find(d == max(d)));

    Q = F1 * x(1) + F2 * x(2);

    R(:, :, j) = Q * V;
    for i = 1:K
      c(i, j) = 1/2 * trace(R(:, :, j) * M(:, :, i)');
    end

    % Measure residual.
    y(j) = norm(M_mat - kron(c(:, j)', R(:, :, j)), 'fro') ^ 2;
  end

  if y(2) < y(1)
    j = 2;
  else
    j = 1;
  end

  R = R(:, :, j);
  c = c(:, j);

  M_mat = kron(c', R);
  M = reshape(M_mat, [2, 3, K]);
end
