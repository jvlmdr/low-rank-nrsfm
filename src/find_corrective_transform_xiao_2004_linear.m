% Parameters:
% M_hat -- 2F x 3K
% subset -- Vector of K frame indices to use as bases.
%
% Returns:
% G -- 3K x 3K corrective transform.
% Rs -- 2 x 3 x F rotation matrices.
% C -- F x K coefficient matrix.

function [G, Rs, C] = find_corrective_transform_xiao_2004_linear(M_hat, subset)
  % M_hat is 2F x 3K.
  F = size(M_hat, 1) / 2;
  K = size(M_hat, 2) / 3;
  n = 3 * K * (3 * K + 1) / 2;

  % Make subset of K frames first.
  not_subset = setdiff(1:F, subset);
  order = [subset(:); not_subset(:)];

  % [2F, 3K] -> [2, F, 3K]
  M_hat = reshape(M_hat, [2, F, 3 * K]);
  % Re-order rows of M_hat.
  M_hat = M_hat(:, order, :);
  % [2, F, 3K] -> [2F, 3K]
  M_hat = reshape(M_hat, [2 * F, 3 * K]);

  % Symmetric parametrization of 3K x 3K matrix.
  H = construct_symmetric(3 * K);

  % Find each corrective triple up to rotation.
  G = zeros(3 * K, 3, K);

  for k = 1:K
    % Build linear system.
    [A, c] = rotation_and_basis_constraints(M_hat, k);
    % Solve linear system.
    q = [A' * A, c; c', 0] \ [zeros(n, 1); 1];
    q = q(1:n);
    % Construct symmetric matrix.
    Q = reshape(H * q, [3 * K, 3 * K]);

    % Extract 3K x 3 matrix by decomposition.
    [V, D] = eig(Q);
    d = diag(D);
    % All eigenvalues should be real because Q is symmetric.
    assert(all(isreal(d)));
    % Take three biggest eigenvalues.
    [d, ind] = sort(d, 'descend');
    V = V(:, ind);
    % Check that Q is rank 3 and positive definite.
    assert(all(d(1:3) > 0));
    % Take decomposition.
    G(:, :, k) = V(:, 1:3) * diag(sqrt(d(1:3)));
  end

  % Determine alignment and signs of coefficients.
  [G, Rs, C] = align_corrective_triples(M_hat, G);

  % Return to original ordering.
  disorder = zeros(F, 1);
  disorder(order) = 1:F;
  Rs = Rs(:, :, disorder);
  C = C(disorder, :);

  % [3K, 3, K] -> [3K, 3K]
  G = reshape(G, [3 * K, 3 * K]);
end
