% [rotations, coeff, M_hat, B_hat, G] = find_rotations_trace(projections, K,
%     lambda)
%
% Parameters:
% projections -- 2 x P x F
% K -- Basis size
% lambda -- Regularization of projections versus trace norm
%
% Returns:
% rotations -- 2 x 3 x F
% coeff -- F x 1
% M_hat -- 2F x 3K
% B_hat -- 3K x P
% G -- Corrective triple. 3K x 3

function [rotations, coeff, M_hat, B_hat, G] = find_rotations_trace(...
    projections, K, lambda)
  % Find Q = G G', where G is the 3K x 3 matrix which best corrects M_hat.

  [M_hat, B_hat, W] = factorize_projections(projections, K);

  H = construct_symmetric(3 * K);
  n = (3 * K) * (3 * K + 1) / 2;

  % Get constraints that make A orthogonal.
  [A, c] = rotation_constraints(M_hat);
  m = size(A, 2);
  % c' x = d.
  C = c';
  d = 1;

  % Find an estimate using trace-norm minimization.
  Q = find_corrective_gram_matrix(A, c, lambda);
  g = H \ Q(:);
  fprintf('rank(Q) = %d\n', rank(Q));
  fprintf('|| A g || / || g || = %g\n', norm(A * g) / norm(g));

  % Project down on to rank-3.
  Q = project_rank(Q, 3);
  g = H \ Q(:);
  fprintf('rank(Q) = %d\n', rank(Q));
  fprintf('|| A g || / || g || = %g\n', norm(A * g) / norm(g));

  % Non-convex refinement with explicit incorporation of rank constraints.
  Q = refine_corrective_gram_matrix(A, c, Q);
  g = H \ Q(:);
  fprintf('rank(Q) = %d\n', rank(Q));
  fprintf('|| A g || / || g || = %g\n', norm(A * g) / norm(g));

  % Extract corrective transform from Gram matrix.
  [V, D] = eig(Q);
  d = diag(D);
  [d, order] = sort(d, 'descend');
  V = V(:, order);
  G = V(:, 1:3) * diag(sqrt(d(1:3)));

  % Apply corrective transform.
  M = M_hat * G;

  % Find nearest rotation matrices and coefficients.
  [rotations, coeff] = nearest_scaled_rotation_matrices(M);
end
