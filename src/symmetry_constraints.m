function A = symmetry_constraints(n)

% X = X^T is equivalent to vec(X) - vec(X^T) = A vec(X) = 0.
A = eye(n * n);
order = reshape(1:(n * n), [n, n])';
A = A(order, :);
A = eye(n * n) - A;

% Take linearly independent rows only.
subset = reshape(1:(n * n), [n, n]);
subset = tril(subset, -1);
subset = squareform(subset);
A = A(subset, :);

end
