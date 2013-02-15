function S_hat = align_structure(S_hat, S)
  F = size(S, 1) / 3;
  P = size(S, 2);

  % [3F, P] -> [3, F, P] -> [3, FP] -> [FP, 3]
  A = reshape(S_hat, [3, F * P])';
  B = reshape(S, [3, F * P])';

  R = procrustes(A, B);
  A = A * R;

  % [FP, 3] -> [3, FP] -> [3, F, P] -> [3F, P]
  S_hat = reshape(A', [3 * F, P]);
end
