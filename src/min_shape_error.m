% Parameters:
% S, X -- 3 x P x F non-rigid shape
% S is the reference, to which X is aligned.

function residual = min_shape_error(S, X)
  assert(ndims(S) == 3);
  assert(size(S, 1) == 3);
  P = size(S, 2);
  F = size(S, 3);
  assert(ndims(X) == 3);
  assert(size(X, 2) == P);
  assert(size(X, 3) == F);

  X = align_shapes_similarity(S, X);

  residual = 0;
  total = 0;

  for t = 1:F
    S_t = S(:, :, t);
    X_t = X(:, :, t);
    D_t = S_t - X_t;
    residual = residual + mean(sqrt(sum(D_t .* D_t)));
    total = total + mean(sqrt(sum(S_t .* S_t)));
  end

  residual = residual / total;
end
