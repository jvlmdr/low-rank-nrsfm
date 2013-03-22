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
    B = inv(G) * B_hat;
    S = kron(C, eye(3)) * B;

    % Compute reprojection error.
    err1 = norm(W - R * S, 'fro') / norm(W, 'fro');

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

    for t = 1:F
      W_t = projections(:, :, t);
      S_t = structure(:, :, t);

      % minimize || S_t' Q - W_t' || s.t. Q' Q = I
      Q = procrustes(S_t', W_t');
      err_new = norm(S_t' * Q - W_t', 'fro');

      Q_old = rotations(:, :, t)';
      err_old = norm(S_t' * Q_old - W_t', 'fro');

      % Protect against numerical error.
      if err_new > err_old
        Q = Q_old;
      end

      rotations(:, :, t) = Q';
    end

    % Compute reprojection error.
    R = block_diagonal_cameras(rotations);
    err2 = norm(W - R * S, 'fro') / norm(W, 'fro');

    fprintf('%6d: %12g %12g\n', i, err1, err2);
  end
end
