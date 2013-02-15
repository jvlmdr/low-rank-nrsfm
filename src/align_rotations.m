function R_hat = align_rotations(R_hat, R)
  G = procrustes(R_hat, R);
  R_hat = R_hat * G;
end
