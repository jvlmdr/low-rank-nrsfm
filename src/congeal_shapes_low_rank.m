% [structure, rotations] = congeal_shapes_low_rank(structure, K, tol, max_iter)
%
% Parameters:
% structure -- 3 x P x F. Assume zero-mean.
% tol -- Minimum relative change to keep iterating.
% max_iter -- Maximum number of iterations.
%
% Returns:
% structure -- 3 x P x F
% rotations -- 3 x 3 x F

function [structure, rotations] = congeal_shapes_low_rank(structure, K, tol, ...
    max_iter)
  % Take a copy.
  init = structure;

  P = size(structure, 2);
  F = size(structure, 3);
  num_iter = 0;
  converged = false;

  while ~converged && num_iter < max_iter
    prev_structure = structure;

    for i = 1:F
      % Align all other shapes to this shape.
      A = repmat(structure(:, :, i), [1, 1, F - 1]);
      B = structure(:, :, [1:(i - 1), (i + 1):F]);
      B = k_reshape(structure_to_matrix(B), 3);
      B = project_rank(B, K);
      B = structure_from_matrix(k_unreshape(B, 3));

      A = A(:, :);
      B = B(:, :);
      R = procrustes(A', B')';

      structure(:, :, i) = R * structure(:, :, i);
    end

    % Check difference.
    % Normalize back to same global transform.
    A = structure(:, :);
    B = prev_structure(:, :);
    R = procrustes(A', B')';
    A = R * A;
    structure = reshape(A, size(structure));
    difference = norm(A(:) - B(:)) / norm(B(:));

    if difference < tol
      converged = true;
    end

    fprintf('%4d: %8g\n', num_iter, difference);
    num_iter = num_iter + 1;
  end

  % Calculate rotations retro-actively.
  rotations = zeros(3, 3, F);
  for t = 1:F
    A = init(:, :, t);
    B = structure(:, :, t);
    R = procrustes(A', B')';
    rotations(:, :, t) = R;
  end
end
