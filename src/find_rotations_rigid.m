function rotations = find_rotations_rigid(projections)
  F = size(projections, 3);
  P = size(projections, 2);

  W = projections_to_matrix(projections);
  [R_hat, S_hat] = factorize_projections(projections, 1);

  [A, c] = rotation_constraints(R_hat);
  C = c';
  d = 1;

  % min x' P x  s.t.  A x = b, X >= 0

  [~, subset] = construct_symmetric(3);

  cvx_begin sdp;
    variable Q(3, 3) symmetric;
    q = Q(subset)';
    
    minimize sum_square(A * q)
    subject to
      % Q = G G', so it must be positive semidefinite.
      Q >= 0;
      % ||i_1|| == 1.
      C * q - d == 0;
  cvx_end;

  % Extract corrective transform from Gram matrix.
  [V, D] = eig(Q);
  d = diag(D);
  % Sort largest to smallest.
  [d, order] = sort(d, 'descend');
  V = V(:, order);
  % Should be positive. Ensure.
  d = max(d, 0);
  % Q = G G'
  G = V * diag(sqrt(d));

  % Apply corrective transform.
  R = R_hat * G;
  S = inv(G) * S_hat;

  [rotations, scales] = nearest_scaled_rotation_matrices(R);
end
