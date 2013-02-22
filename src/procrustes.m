% Returns orthonormal R which minimizes || A R - B ||_F.
% http://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem

function R = procrustes(A, B)
  M = A' * B;
  [U, S, V] = svd(M);
  R = U * V';
end
