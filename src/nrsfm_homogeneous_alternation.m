% [structure, rotations, basis, coeff] =
%   nrsfm_homogeneous_alternation(projections, rotations, K)
%
% Parameters:
% projections -- 2 x P x F
% rotations -- 2 x 3 x F
% basis -- 3 x K x P
% K -- Shape basis size
%
% Returns:
% structure -- 3 x P x F
% rotations -- 2 x 3 x F
% basis -- 3 x K x P
% coeff -- K x F

function [structure, rotations, basis, coeff] = ...
    nrsfm_homogeneous_alternation(projections, rotations, basis, num_iter)
  P = size(projections, 2);
  F = size(projections, 3);
  K = size(basis, 2);

  [M_hat, B_hat, W] = factorize_projections(projections, K);
  R = block_diagonal_cameras(rotations);

  for i = 1:num_iter
    % [3, K, P] -> [3K, P]
    B = basis;
    B = reshape(B, [3 * K, P]);

    [G, C] = find_corrective_matrix(M_hat, R, B);
    M = M_hat * G;
    W_hat = M * B;
    err1 = norm(W_hat - R * kron(C, eye(3)) * B, 'fro') / norm(W_hat, 'fro');

    B = inv(G) * B_hat;
    S = kron(C, eye(3)) * B;

    % Compute reprojection error.
    %err1 = norm(W - R * S, 'fro') / norm(W, 'fro');

    % [3K, P] -> [3, K, P]
    basis = B;
    basis = reshape(basis, [3, K, P]);
    % [F, K] -> [K, F]
    coeff = C';
    % Normalize the K basis vectors.
    [basis, coeff] = normalize_basis(basis, coeff);

    % [3F, P] -> [3, F, P] -> [3, P, F]
    structure = S;
    structure = reshape(structure, [3, F, P]);
    structure = permute(structure, [1, 3, 2]);
    % [2F, P] -> [2, F, P] -> [2, P, F]
    reprojections = W_hat;
    reprojections = reshape(reprojections, [2, F, P]);
    reprojections = permute(reprojections, [1, 3, 2]);

    %rotations = find_non_worse_cameras(projections, structure, rotations);
    rotations = find_non_worse_cameras(reprojections, structure, rotations);

    % Compute reprojection error.
    R = block_diagonal_cameras(rotations);
    %err2 = norm(W - R * S, 'fro') / norm(W, 'fro');
    err2 = norm(M - R * kron(C, eye(3)), 'fro') / norm(M, 'fro');

    fprintf('%6d: %12g %12g\n', i, err1, err2);
  end
end
