function R = procrustes(A, B)
  M = A' * B;
  [U, S, V] = svd(M);
  R = U * V';
end
