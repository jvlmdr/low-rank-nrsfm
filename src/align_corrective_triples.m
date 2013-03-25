% Parameters:
% M_hat -- 2F x 3K, obtained from SVD.
% G -- 3K x 3 x K, corrective triple for each basis.
%
% Returns:
% G -- 3K x 3 x K
% rotations -- 2 x 3 x F
% coeff -- K x F

function [G, rotations, coeff] = align_corrective_triples(M_hat, G)
  % M_hat is 2F x 3K.
  F = size(M_hat, 1) / 2;
  K = size(M_hat, 2) / 3;

  % [3K, 3, K] -> [3K, 3K]
  G = reshape(G, [3 * K, 3 * K]);

  % Apply all unaligned triples to the uncorrected M.
  M = M_hat * G;

  % [3K, 3K] -> [3K, 3, K]
  G = reshape(G, [3 * K, 3, K]);
  % [2F, 3K] -> [2, F, 3, K] -> [F, K, 2, 3]
  M = reshape(M, [2, F, 3, K]);
  M = permute(M, [2, 4, 1, 3]);

  % Extract coefficients (up to sign) from M.
  % Find scale such that the mean of norm(i) and norm(j) = 1.
  C = mean(sqrt(sum(M .^ 2, 4)), 3);

  % [F, K, 2, 3] -> [2, 3, F, K]
  M = permute(M, [3, 4, 1, 2]);

  for k = 2:K
    % Pick frame to align.
    cc = C(:, 1) .* C(:, k);
    [max_cc, u] = max(cc);

    % Align camera u from triple k to that of triple 1.
    % Assume C(u, 1) > 0 and C(u, k) > 0.
    % Scale doesn't matter to orthogonal Procrustes since it's simply an SVD.
    M_u1 = M(:, :, u, 1);
    M_uk = M(:, :, u, k);
    [E(:, :, 1), E(:, :, 2)] = ambiguous_procrustes(M_uk, M_u1);

    % Evaluate both transforms.
    residuals = zeros(F, 2);
    signs = zeros(F, 2);

    for j = 1:2
      % Take a copy.
      R_1 = M(:, :, :, 1);
      R_k = M(:, :, :, k);

      % [2, 3, F] -> [2, F, 3] -> [2F, 3]
      R_k = permute(R_k, [1, 3, 2]);
      R_k = reshape(R_k, [2 * F, 3]);

      % Apply transform from Procrustes to all cameras.
      R_k = R_k * E(:, :, j);

      % [2F, 3] -> [2, F, 3] -> [2, 3, F]
      R_k = reshape(R_k, [2, F, 3]);
      R_k = permute(R_k, [1, 3, 2]);

      % Scale each rotation by the other basis' scale to avoid division.
      R_1 = bsxfun(@times, R_1, reshape(C(:, k), [1, 1, F]));
      R_k = bsxfun(@times, R_k, reshape(C(:, 1), [1, 1, F]));

      % Determine sign of each C(t, k).
      for t = 1:F
        positive_residual = norm(R_1(:, :, t) - R_k(:, :, t));
        negative_residual = norm(R_1(:, :, t) + R_k(:, :, t));

        if negative_residual < positive_residual
          signs(t, j) = -1;
        else
          signs(t, j) = 1;
        end

        residuals(t, j) = min(positive_residual, negative_residual);
      end
    end

    % Pick between ambiguous Procrustes solutions.
    residuals = mean(residuals);
    fprintf('Ambiguous Procrustes residuals: %g, %g\n', min(residuals), ...
        max(residuals));
    [residual, arg] = min(residuals);
    E = E(:, :, arg);
    signs = signs(:, arg);

    C(:, k) = signs .* C(:, k);

    % Re-extract cameras.
    R_1 = M(:, :, :, 1);
    R_k = M(:, :, :, k);

    % Scale each rotation by the other basis' scale to avoid division.
    R_1 = bsxfun(@times, R_1, reshape(C(:, k), [1, 1, F]));
    R_k = bsxfun(@times, R_k, reshape(C(:, 1), [1, 1, F]));

    % [2, 3, F] -> [2, F, 3] -> [2F, 3]
    R_1 = permute(R_1, [1, 3, 2]);
    R_1 = reshape(R_1, [2 * F, 3]);
    R_k = permute(R_k, [1, 3, 2]);
    R_k = reshape(R_k, [2 * F, 3]);

    % After fixing signs, align all cameras.
    E = procrustes(R_k, R_1);

    % [2, 3, F, K] -> [2, F, 3, K] -> [2F, 3, K]
    M = permute(M, [1, 3, 2, 4]);
    M = reshape(M, [2 * F, 3, K]);

    % Apply transform to triple.
    G(:, :, k) = G(:, :, k) * E;
    M(:, :, k) = M(:, :, k) * E;

    % [2F, 3, K] -> [2, F, 3, K] -> [2, 3, F, K]
    M = reshape(M, [2, F, 3, K]);
    M = permute(M, [1, 3, 2, 4]);
  end

  % [2, 3, F, K] -> [2, K, 3, F] -> [2K, 3, F]
  M = permute(M, [1, 4, 2, 3]);
  M = reshape(M, [2 * K, 3, F]);

  rotations = zeros(2, 3, F);

  % Now extract rotations from the [c_t1 R_t ... c_tK R_t] rows.
  for t = 1:F
    % Align scaled identity matrices to scaled rotations.
    A = kron(C(t, :)', eye(2, 3));
    R_t = procrustes(A, M(:, :, t));
    rotations(:, :, t) = R_t(1:2, :);
  end

  % [F, K] -> [K, F]
  coeff = C';
end
