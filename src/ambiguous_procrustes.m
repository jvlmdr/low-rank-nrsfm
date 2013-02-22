% Returns orthonormal R which minimizes || A R - B ||_F.
% http://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem
%
% When the matrix M = A' * B is rank 2, the last vector is up to sign.

function [R1, R2] = ambiguous_procrustes(A, B)
  M = A' * B;
  [U, S, V] = svd(M);
  R1 = U * V';
  R2 = U * diag([1, 1, -1]) * V';
end
