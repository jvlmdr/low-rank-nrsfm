function [R, k] = nearest_scaled_rotation_matrix(P)
  [U, S, V] = svd(P, 'econ');
  R = U * V';
  k = R(:) \ P(:);
end
