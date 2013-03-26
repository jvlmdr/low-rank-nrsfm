% Parameters:
% S, X -- 3 x P x F non-rigid shape
% S is the reference, to which X is aligned.

function residual = min_total_shape_error(S, X)
  assert(ndims(S) == 3);
  assert(size(S, 1) == 3);
  P = size(S, 2);
  F = size(S, 3);
  assert(ndims(X) == 3);
  assert(size(X, 2) == P);
  assert(size(X, 3) == F);

  residuals = zeros(F, 1);
  totals = zeros(F, 1);

  for t = 1:F
    [residuals(t), totals(t)] = min_shape_error(S, X);
  end

  residual = sum(residual) / sum(total);
end
