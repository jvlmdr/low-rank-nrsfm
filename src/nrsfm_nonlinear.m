% Parameters:
% projections -- 2 x P x F
% rotations -- 2 x 3 x F
% basis -- 3 x K x P
% coeff -- K x F
%
% Returns:
% structure -- 3 x P x F
% rotations -- 2 x 3 x F
% basis -- 3 x K x P
% coeff -- K x F

function [structure, rotations, basis, coeff] = nrsfm_nonlinear(projections, ...
    rotations, basis, coeff, max_iter, tol)
  assert(ndims(projections) == 3);
  assert(size(projections, 1) == 2);
  P = size(projections, 2);
  F = size(projections, 3);
  K = size(basis, 2);

  assert(ndims(rotations) == 3);
  assert(all(size(rotations) == [2, 3, F]));

  assert(ndims(basis) == 3);
  assert(all(size(basis) == [3, K, P]));

  assert(ndims(coeff) == 2);
  assert(all(size(coeff) == [K, F]));

  % Convert rotations to quaternions.
  quaternions = zeros(4, F);
  for t = 1:F
    R_t = rotations(:, :, t);
    R_t = [R_t; cross(R_t(1, :), R_t(2, :))];
    quaternions(:, t) = rot2quat(R_t);
  end

  % Compute initial residual for debug purposes. Should reassure us that
  % quaternion conversion worked and matrices are being addressed correctly.
  W = projections_to_matrix(projections);
  R = block_diagonal_cameras(rotations);
  structure = compose_structure(basis, coeff);
  S = structure_to_matrix(structure);
  fprintf('Initial residual: %g\n', 1/2 * norm(W - R * S, 'fro') ^ 2);

  % Solve.
  [quaternions, basis, coeff] = nrsfm_nonlinear_mex(projections, ...
      quaternions, basis, coeff, max_iter, tol);

  structure = compose_structure(basis, coeff);

  % Convert back to rotation matrices.
  rotations = zeros(2, 3, F);
  for t = 1:F
    R_t = quat2rot(quaternions(:, t));
    rotations(:, :, t) = R_t(1:2, :);
  end

  % Compute final residual for debug purposes.
  R = block_diagonal_cameras(rotations);
  S = structure_to_matrix(structure);
  fprintf('Final residual: %g\n', 1/2 * norm(W - R * S, 'fro') ^ 2);
end
