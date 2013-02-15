function R = find_rotations(P_hat, lambda)
  % Find G = Q Q', where Q is the 3K x 3 matrix which best corrects P_hat.

  F = size(P_hat, 1) / 2;
  K = size(P_hat, 2) / 3;

  H = construct_symmetric(3 * K);
  n = (3 * K) * (3 * K + 1) / 2;

  % Get constraints that make A orthogonal.
  [A, c] = rotation_constraints(P_hat);
  m = size(A, 2);
  % c' x = d.
  C = c';
  d = 1;

  % Find an estimate using trace-norm minimization.
  G = find_corrective_gram_matrix(A, c, lambda);
  g = H \ G(:);
  fprintf('rank(G) = %d\n', rank(G));
  fprintf('|| A g || / || g || = %g\n', norm(A * g) / norm(g));

  % Project down on to rank-3.
  G = project_rank(G, 3);
  g = H \ G(:);
  fprintf('rank(G) = %d\n', rank(G));
  fprintf('|| A g || / || g || = %g\n', norm(A * g) / norm(g));

  % Non-convex refinement with explicit incorporation of rank constraints.
  G = refine_corrective_gram_matrix(A, c, G);
  g = H \ G(:);
  fprintf('rank(G) = %d\n', rank(G));
  fprintf('|| A g || / || g || = %g\n', norm(A * g) / norm(g));

  % Extract corrective transform from Gram matrix.
  [V, D] = eig(G);
  d = diag(D);
  [d, order] = sort(d, 'descend');
  V = V(:, order);
  Q = V(:, 1:3) * diag(sqrt(d(1:3)));

  % Apply corrective transform.
  P = P_hat * Q;

  % Find nearest rotation matrices and coefficients.
  [R, c] = nearest_scaled_rotation_matrices(P);

%  % Negate where required.
%  for t = 1:F - 1
%    R1 = R(:, :, t);
%    R2 = R(:, :, t + 1);
%    k1 = cross(R1(1, :), R1(2, :));
%    k2 = cross(R2(1, :), R2(2, :));
%    R1 = [R1; k1];
%    R2 = [R2; k2];
%
%    theta(t) = acos((trace(R2' * R1) - 1) / 2) * 180 / pi;
%
%    if theta(t) > 90
%      R(:, :, t + 1) = -R(:, :, t + 1);
%    end
%  end

  % Negate where required.
  for t = 1:F - 1
    k1 = cross(R(1, :, t), R(2, :, t));
    k2 = cross(R(1, :, t + 1), R(2, :, t + 1));

    if dot(k1, k2) < 0
      warning('Flipping sign of cameras');
      R(:, :, t + 1) = -R(:, :, t + 1);
    end
  end

  % [2, 3, F] -> [2, F, 3] -> [2F, 3]
  R = reshape(permute(R, [1, 3, 2]), [2 * F, 3]);
end
