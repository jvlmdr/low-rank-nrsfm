% Parameters:
% S, X -- 3 x P shape
%
% S is the reference, to which X is aligned.

function X = align_shapes_similarity(S, X)
  assert(ndims(S) == 2);
  assert(size(S, 1) == 3);
  P = size(S, 2);
  assert(ndims(X) == 2);
  assert(all(size(X) == size(S)));

  % [3, P] -> [P, 3]
  S = S';
  X = X';

  % Remove centroid.
  centroid_X = mean(X);
  centroid_S = mean(S);
  X = bsxfun(@minus, X, centroid_X);
  S = bsxfun(@minus, S, centroid_S);

  % Normalize scale.
  x = sqrt(mean(sum(X .^ 2, 2), 1));
  s = sqrt(mean(sum(S .^ 2, 2), 1));
  if x == 0
    s = 0;
  else
    s = s / x;
  end
  X = s * X;

  % Align each frame individually.
  R = procrustes(X, S);
  X = X * R;

  % Restore centroid.
  X = bsxfun(@plus, X, centroid_S);

  % [P, 3] -> [3, P]
  X = X';
end
