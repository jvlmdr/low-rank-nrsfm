% Returns the rotation R which minimizes || A R - B ||_F.
% http://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem

function R = rotation_procrustes(A, B)
  M = A' * B;
  [U, S, V] = svd(M);
  R = U * diag([1, 1, sign(det(U * V'))]) * V';
end
