% [structure, rotations, basis, coeff] =
%   nrsfm_nullspace_alternation(projections, rotations, K, max_iter)
%
% Parameters:
% projections -- 2 x P x F
% rotations -- 2 x 3 x F
% K -- Shape basis size
% max_iter -- Number of iterations of alternation to do
%
% Returns:
% structure -- 3 x P x F
% rotations -- 2 x 3 x F
% basis -- 3 x K x P
% coeff -- K x F

function [structure, rotations, basis, coeff] = nrsfm_nullspace_alternation(...
    projections, rotations, K, max_iter)
  P = size(projections, 2);
  F = size(projections, 3);

  [M_hat, B_hat, W] = factorize_projections(projections, K);

  for i = 1:max_iter
    R = block_diagonal_cameras(rotations);
    [G, C] = find_corrective_matrix_nullspace(M_hat, R);
    B = inv(G) * B_hat;
    S = kron(C, eye(3)) * B;

    err1 = norm(W - R * S, 'fro') / norm(W, 'fro');

    % [3F, P] -> [3, F, P] -> [3, P, F]
    structure = S;
    structure = reshape(structure, [3, F, P]);
    structure = permute(structure, [1, 3, 2]);

    rotations = find_non_worse_cameras(projections, structure, rotations);

    R = block_diagonal_cameras(rotations);
    err2 = norm(W - R * S, 'fro') / norm(W, 'fro');

    fprintf('%6d: %12g %12g\n', i, err1, err2);
  end

  % [3K, P] -> [3, K, P]
  basis = B;
  basis = reshape(basis, [3, K, P]);
  % [F, K] -> [K, F]
  coeff = C';
end
