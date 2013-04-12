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
    % Evaluate both transforms.
    residuals = zeros(2, F);
    signs = zeros(F, 2, F);

    % Align using frame u.
    for u = 1:F
      cc = C(:, 1) .* C(:, k);

      % Align camera u from triple k to that of triple 1.
      % Assume C(u, 1) > 0 and C(u, k) > 0.
      % Scale doesn't matter to orthogonal Procrustes since it's simply an SVD.
      M_u1 = M(:, :, u, 1);
      M_uk = M(:, :, u, k);
      [E(:, :, 1, u), E(:, :, 2, u)] = ambiguous_procrustes(M_uk, M_u1);

      for j = 1:2
        % Take a copy.
        R_1 = M(:, :, :, 1);
        R_k = M(:, :, :, k);

        % [2, 3, F] -> [2, F, 3] -> [2F, 3]
        R_k = stack_cameras(R_k);
        % Apply transform from Procrustes to all cameras.
        R_k = R_k * E(:, :, j, u);
        % [2F, 3] -> [2, F, 3] -> [2, 3, F]
        R_k = unstack_cameras(R_k);

        % Scale each rotation by the other basis' scale to avoid division.
        R_1 = bsxfun(@times, R_1, reshape(C(:, k), [1, 1, F]));
        R_k = bsxfun(@times, R_k, reshape(C(:, 1), [1, 1, F]));

        % Determine sign of C(t, k) for each t.
        for t = 1:F
          positive_residual = norm(R_1(:, :, t) - R_k(:, :, t), 'fro');
          negative_residual = norm(R_1(:, :, t) + R_k(:, :, t), 'fro');

          if negative_residual < positive_residual
            signs(t, j, u) = -1;
          else
            signs(t, j, u) = 1;
          end
        end

        % Apply signs to coefficients.
        c_k = signs(:, j, u) .* C(:, k);
        % Re-extract cameras.
        R_1 = M(:, :, :, 1);
        R_k = M(:, :, :, k);

        % Scale each rotation by the other basis' scale to avoid division.
        R_1 = bsxfun(@times, R_1, reshape(c_k, [1, 1, F]));
        R_k = bsxfun(@times, R_k, reshape(C(:, 1), [1, 1, F]));

        % [2, 3, F] -> [2, F, 3] -> [2F, 3]
        R_1 = stack_cameras(R_1);
        R_k = stack_cameras(R_k);

        % After fixing signs, align all cameras.
        E(:, :, j, u) = procrustes(R_k, R_1);

        % Compute residual.
        residuals(j, u) = norm(R_k * E(:, :, j, u) - R_1, 'fro');
      end
    end

    % Find best transform.
    [min_residual, index] = min(residuals(:));
    [j_star, u_star] = ind2sub(size(residuals), index);
    E = E(:, :, index);
    signs = signs(:, index);

    % Apply best signs.
    C(:, k) = C(:, k) .* signs;

    % Apply transform to triple.
    G(:, :, k) = G(:, :, k) * E;

    % [2, 3, F, K] -> [2, F, 3, K] -> [2F, 3, K]
    M = permute(M, [1, 3, 2, 4]);
    M = reshape(M, [2 * F, 3, K]);
    % Apply transform to motion matrix.
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
