% Constructs a vectorized symmetric matrix from the unique entries.
% Maps n * (n + 1) / 2 -> n * n.
% reshape(construct_symmetric(n) * x, [n, n]) gives the symmetric matrix.

function [A, subset] = construct_symmetric(n)
  A1 = speye(n * n);

  order = reshape(1:(n * n), [n, n])';
  A2 = A1(order, :);

  A = A1 + A2;

  % Prevent diagonal elements from being doubled.
  A = double(A ~= 0);

  % Take linearly independent rows only.
  subset = reshape(1:(n * n), [n, n])';
  subset = tril(subset, -1);
  subset = squareform(subset);
  % Invert subset.
  subset = setdiff(1:(n * n), subset);

  A = A(:, subset);
end
