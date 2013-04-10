% [structure, basis, coeff] = find_structure_approx_low_rank_nonlinear(
%     projections, rotations, basis, coeff, max_iter, tol)
%
% Parameters:
% projections -- 2 x P x F
% structure -- 3 x P x F
% basis -- 3 x K x P
% coeff -- K x F
%
% Returns:
% structure -- 3 x P x F
% rotations -- 2 x 3 x F
% basis -- 3 x K x P
% coeff -- K x F

function structure = find_structure_approx_low_rank_nonlinear(projections, ...
    structure, rotations, basis, coeff, lambda, max_iter, tol)
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
  d = structure - compose_structure(basis, coeff);
  r = 1/2 * projection_error(projections, structure, rotations) ^ 2 + ...
      lambda / 2 * norm(d(:)) ^ 2;
  fprintf('Initial residual: %g\n', r);
  r_init = r;

  % Solve.
  if exist('find_structure_approx_low_rank_nonlinear_mex', 'file')
    [structure, basis, coeff] = find_structure_approx_low_rank_nonlinear_mex(...
        projections, structure, quaternions, basis, coeff, lambda, max_iter, ...
        tol);
  else
    error('Not compiled');
  end

  % Convert back to rotation matrices.
  rotations = zeros(2, 3, F);
  for t = 1:F
    R_t = quat2rot(quaternions(:, t));
    rotations(:, :, t) = R_t(1:2, :);
  end

  % Compute final residual for debug purposes.
  d = structure - compose_structure(basis, coeff);
  r = 1/2 * projection_error(projections, structure, rotations) ^ 2 + ...
      lambda / 2 * norm(d(:)) ^ 2;
  fprintf('Initial residual: %g\n', r_init);
  fprintf('Final residual: %g\n', r);
  r_init = r;
end
