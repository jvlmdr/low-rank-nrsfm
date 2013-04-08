% Returns orthonormal R which minimizes || A R - B ||_F.
% http://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem

function R = procrustes(A, B)
  assert(all(size(A) == size(B)));
  M = A' * B;
  [U, S, V] = svd(M);
  R = U(:, 1:min(size(S))) * V';
end
