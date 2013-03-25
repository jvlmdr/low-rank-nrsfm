% Solves
% arg min_{M, B} ||W - MB||_F  s.t.  M \in {motion matrices}
%
% [structure, rotations, basis, coeff] =
%   nrsfm_metric_projections(projections, rotations, basis, max_iter)
%
% Parameters:
% projections -- 2 x P x F
% rotations -- 2 x 3 x F
% coeff -- K x F
% max_iter -- Number of iterations of alternation to do
%
% Returns:
% structure -- 3 x P x F
% rotations -- 2 x 3 x F
% basis -- 3 x K x P
% coeff -- K x F

function [structure, rotations, basis, coeff] = nrsfm_metric_projections(...
    projections, rotations, coeff, max_iter)
  % Introduce the auxiliary variable N.
  % arg min_{M, N, B} ||W - MB||_F s.t. N \in {motion matrices}, M - N = 0

  P = size(projections, 2);
  F = size(projections, 3);
  K = size(coeff, 1);

  % Compose motion matrix.
  % [K, F] -> [F, K]
  C = coeff';
  % [2, 3, F] -> [2F, 3F]
  R = block_diagonal_cameras(rotations);
  M = R * kron(C, eye(3));
  % [2F, 3K] -> [2, F, 3, K]
  M = reshape(M, [2, F, 3, K]);

  % [2, P, F] -> [2, F, P] -> [2F, P]
  W = reshape(permute(projections, [1, 3, 2]), [2 * F, P]);

  B = zeros(3 * K, P);

  for i = 1:max_iter
    % Project on to motion manifold.
    % [2, F, 3, K] -> [2, 3, K, F]
    M = permute(M, [1, 3, 4, 2]);
    for t = 1:F
      [M(:, :, :, t), c_t, R_t] = project_motion_manifold(M(:, :, :, t));
      % For output.
      coeff(:, t) = c_t;
      rotations(:, :, t) = R_t;
    end
    % [2, 3, K, F] -> [2, F, 3, K]
    M = permute(M, [1, 4, 2, 3]);

    % [2, F, 3, K] -> [2F, 3K]
    M = reshape(M, [2 * F, 3 * K]);
    err1 = norm(W - M * B, 'fro');

    % Minimize projection error.
    B = M \ W;
    err2 = norm(W - M * B, 'fro');
    % [3K, P] -> [3, K, P]
    basis = reshape(B, [3, K, P]);

    % Minimize projection error.
    M = W / B;
    err3 = norm(W - M * B, 'fro');
    % [2F, 3K] -> [2, F, 3, K]
    M = reshape(M, [2, F, 3, K]);

    fprintf('%6d: %8g %8g %8g\n', i, err1, err2, err3);
    %fprintf('%6d\n', i);
  end

  structure = compose_structure(basis, coeff);
end
