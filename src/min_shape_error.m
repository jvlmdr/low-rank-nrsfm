% Parameters:
% S, X -- 3 x P shape
%
% S is the reference, to which X is aligned.

function [residual, total] = min_shape_error(S, X)
  assert(ndims(S) == 2);
  assert(size(S, 1) == 3);
  P = size(S, 2);
  assert(ndims(X) == 2);
  assert(all(size(X) == size(S)));

  X = align_shapes_similarity(S, X);

  D = S - X;
  residual = mean(sqrt(sum(D .^ 2)));
  total = mean(sqrt(sum(S .^ 2)));
end
