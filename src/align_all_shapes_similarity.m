% Parameters:
% S, X -- 3 x P x F non-rigid shape
%
% S is the reference, to which X is aligned.

function X = align_all_shapes_similarity(S, X)
  assert(ndims(S) == 3);
  assert(size(S, 1) == 3);
  P = size(S, 2);
  F = size(S, 3);
  assert(ndims(X) == 3);
  assert(size(X, 2) == P);
  assert(size(X, 3) == F);

  for t = 1:F
    X(:, :, t) = align_shapes_similarity(S(:, :, t), X(:, :, t));
  end
end
