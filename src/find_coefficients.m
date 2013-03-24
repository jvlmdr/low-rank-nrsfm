% coeff = find_coefficients(projections, rotations, basis)
%
% Parameters:
% projections -- 2 x P x F
% rotations -- 2 x 3 x F
% basis -- 3 x K x P
%
% Returns:
% coeff -- K x F

function coeff = find_coefficients(projections, rotations, basis)
  P = size(projections, 2);
  F = size(projections, 3);
  K = size(basis, 2);

  coeff = zeros(K, F);

  % Get 3P x K basis matrix.
  % [3, K, P] -> [3, P, K] -> [3P, K]
  B = reshape(permute(basis, [1, 3, 2]), [3 * P, K]);
  % Actually make that a 3 x PK basis matrix.
  B = reshape(B, [3, P * K]);

  % Solve w_t = kron(I, R_t) B c_t for each t.
  for t = 1:F
    R_t = rotations(:, :, t);

    % More efficient to compute R_t * reshape(B) than kron(I, R_t) * B.
    A_t = R_t * B;
    % [2, PK] -> [2P, K]
    A_t = reshape(A_t, [2 * P, K]);

    w_t = projections(:, :, t);
    w_t = w_t(:);

    c_t = A_t \ w_t;
    coeff(:, t) = c_t;
  end
end
