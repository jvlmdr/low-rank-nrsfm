function rotations = find_rotations_rigid(projections)
  F = size(projections, 3);
  P = size(projections, 2);

  W = projections_to_matrix(projections);
  [R_hat, S_hat] = factorize_projections(projections, 1);

  [A, C] = rotation_constraints(R_hat);
  %C = C(1, :);
  d = 1;

  [H, subset] = construct_symmetric(3);

%  % First solve using Lagrange multipliers, sans semidefinite constraint.
%  x = [A' * A, C(1, :)'; C(1, :), 0] \ [zeros(6, 1); d];
%  % Extract primal variables only.
%  x = x(1:6);
%  % Build symmetric matrix.
%  X = reshape(H * x, [3, 3]);
%  % Check if matrix is positive semidefinite.
%  lambda = eig(X);
%  if any(and(lambda < 0, abs(lambda) > 1e-6 * max(abs(lambda))))
%    warning('Unconstrained solution not positive semidefinite');
%    disp(sort(lambda, 'descend'));
%  end

  % Solve using CVX.
  cvx_begin sdp;
    cvx_solver sedumi;

    variable Q(3, 3) symmetric;
    q = Q(subset)';
    
    minimize sum_square(A * q)
    subject to
      % Q = G G', so it must be positive semidefinite.
      Q >= 0;
      % ||i_1|| == 1.
      %C(1, :) * q == d;
      C * q >= d;
  cvx_end;

  if isinf(cvx_optval)
    % Feasibility not satisfied!
    rotations = [];
    return;
  end

  Q = Q / norm(Q(:));
  q = Q(subset)';

  fprintf('rank(Q): %d, |A q|: %g\n', rank(Q), norm(A * q));
  fprintf('Singular values of Gram matrix:\n');
  disp(svd(Q));

%  % Solve using SeDuMi.
%  AA = zeros(size(A, 1), 9);
%  AA(:, subset) = A;
%  b = zeros(size(A, 1), 1);
%  CC = zeros(size(C, 1), 9);
%  CC(:, subset) = C(1, :);
%
%  Q = semidefinite_constrained_least_squares(AA, b, CC, d);

  % Extract corrective transform from Gram matrix.
  [V, D] = eig(Q);
  lambda = diag(D);
  % Factorize Q = G G', zero any small negative values.
  G = V * diag(sqrt(max(lambda, 0)));

  QQ = Q;
  qq = q;
  GG = G;

  % Non-linear refinement.
  fprintf('Non-linear refinement\n');
  row_pairs = projections_from_matrix(R_hat);
  G = refine_corrective_triple_nonlinear(row_pairs, G, 1000, 1e-6);
  Q = G * G';
  % Normalize such that norm(Q(:)) = 1.
  G = G / sqrt(norm(Q(:)));
  Q = G * G';
  q = Q(subset)';
  fprintf('rank(Q): %d, |A q|: %g\n', rank(Q), norm(A * q));
  fprintf('Singular values of Gram matrix:\n');
  disp(svd(Q));

  % Apply corrective transform.
  R = R_hat * G;
  S = inv(G) * S_hat;

  rotations = unstack_cameras(R);
  [rotations, scales] = nearest_scaled_rotation_matrices(rotations);
end
