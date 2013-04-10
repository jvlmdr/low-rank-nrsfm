% [rotations, coeff, M_hat, B_hat, G] = find_rotations_trace(projections, K)
%
% Parameters:
% projections -- 2 x P x F
% K -- Basis size
%
% Returns:
% rotations -- 2 x 3 x F
% coeff -- F x 1
% M_hat -- 2F x 3K
% B_hat -- 3K x P
% G -- Corrective triple. 3K x 3

function [rotations, coeff, M_hat, B_hat, G] = find_rotations_trace(...
    projections, K)
  % Find Q = G G', where G is the 3K x 3 matrix which best corrects M_hat.

  [M_hat, B_hat, W] = factorize_projections(projections, K);

  [~, subset] = construct_symmetric(3 * K);
  n = (3 * K) * (3 * K + 1) / 2;

  % Get constraints that make A orthogonal.
  [A, C] = rotation_constraints(M_hat);
  m = size(A, 2);
  d = 1;

  % Find an estimate using trace-norm minimization.
  fprintf('Solve for Q using trace norm + golden section search\n');
  Q = find_corrective_gram_matrix_golden(A, C, 1e-3);
  q = Q(subset)';
  fprintf('rank(Q): %d, A q / q: %g\n', rank(Q), norm(A * q) / norm(q));
  fprintf('Singular values of Gram matrix\n');
  disp(svd(Q));

  % Project down on to rank-3.
  fprintf('Project on to rank 3\n');
  Q = project_rank(Q, 3);
  q = Q(subset)';
  fprintf('rank(Q): %d, A q / q: %g\n', rank(Q), norm(A * q) / norm(q));

%  % Find an estimate using trace-norm minimization.
%  Q = find_corrective_gram_matrix(A, C, lambda);
%  g = H \ Q(:);
%  fprintf('rank(Q) = %d\n', rank(Q));
%  fprintf('|| A g || / || g || = %g\n', norm(A * g) / norm(g));
%
%  % Project down on to rank-3.
%  Q = project_rank(Q, 3);
%  g = H \ Q(:);
%  fprintf('rank(Q) = %d\n', rank(Q));
%  fprintf('|| A g || / || g || = %g\n', norm(A * g) / norm(g));

%  % Non-convex refinement with explicit incorporation of rank constraints.
%  Q = refine_corrective_gram_matrix(A, C, Q);
%  g = H \ Q(:);
%  fprintf('rank(Q) = %d\n', rank(Q));
%  fprintf('|| A g || / || g || = %g\n', norm(A * g) / norm(g));

  % Extract corrective transform from Gram matrix.
  [V, D] = eig(Q);
  d = diag(D);
  [d, order] = sort(d, 'descend');
  V = V(:, order);
  G = V(:, 1:3) * diag(sqrt(d(1:3)));

  % Non-linear refinement.
  fprintf('Non-linear refinement\n');
  row_pairs = projections_from_matrix(M_hat);
  G = refine_corrective_triple_nonlinear(row_pairs, G, 1000, 1e-6);
  Q = G * G';
  q = Q(subset)';
  fprintf('rank(Q): %d, A q / q: %g\n', rank(Q), norm(A * q) / norm(q));
  fprintf('Singular values of Gram matrix\n');
  disp(svd(Q));

  % Apply corrective transform.
  M = M_hat * G;

  % Find nearest rotation matrices and coefficients.
  [rotations, coeff] = nearest_scaled_rotation_matrices(unstack_cameras(M));
end
