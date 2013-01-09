% diag_vec(n) * x = diag(reshape(x, [n, n])).
% Maps n * n -> n.

function A = diag_vec(n)

A = eye(n * n);
subset = reshape(1:(n * n), [n, n]);
subset = diag(subset);
A = A(subset, :);

end
