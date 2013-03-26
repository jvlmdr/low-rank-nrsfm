% [structure, rotations, basis, coeff] =
%   nrsfm_nullspace_alternation_algebraic(projections, rotations, K, max_iter)
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

function [structure, rotations, basis, coeff] = ...
    nrsfm_nullspace_alternation_algebraic(projections, rotations, K, max_iter)
  P = size(projections, 2);
  F = size(projections, 3);

  [M_hat, B_hat, W] = factorize_projections(projections, K);

  for i = 1:max_iter
    R = block_diagonal_cameras(rotations);
    [G, C] = find_corrective_matrix_nullspace(M_hat, R);

    B = inv(G) * B_hat;
    S = kron(C, eye(3)) * B;
    M = M_hat * G;
    err1 = norm(M - R * kron(C, eye(3)), 'fro');

    % [2F, 3K] -> [2, F, 3K] -> [2, 3K, F]
    target = M;
    target = reshape(target, [2, F, 3 * K]);
    target = permute(target, [1, 3, 2]);

    % [3F, 3K] -> [3, F, 3K] -> [3, 3K, F]
    source = kron(C, eye(3));
    source = reshape(source, [3, F, 3 * K]);
    source = permute(source, [1, 3, 2]);

    % minimize || M_hat_t - R_t kron(c_t', eye(3)) || s.t. Q' Q = I
    rotations = find_non_worse_cameras(target, source, rotations);

    R = block_diagonal_cameras(rotations);
    err2 = norm(M - R * kron(C, eye(3)), 'fro');

    fprintf('%6d: %12g %12g\n', i, err1, err2);
  end

  % [3F, P] -> [3, F, P] -> [3, P, F]
  structure = S;
  structure = reshape(structure, [3, F, P]);
  structure = permute(structure, [1, 3, 2]);

  % [3K, P] -> [3, K, P]
  basis = B;
  basis = reshape(basis, [3, K, P]);
  % [F, K] -> [K, F]
  coeff = C';
end
