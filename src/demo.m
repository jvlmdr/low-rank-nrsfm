clear;
close all;
mocap_file = '../data/akhter-2008/yoga.mat';
K = 8;

load(mocap_file, 'S', 'Rs');
F = size(Rs, 1) / 2;
N = size(S, 2);

% Build block-diagonal rotation matrix.
R = mat2cell(Rs, 2 * ones(F, 1), 3);
R = cellfun(@(X) { sparse(X) }, R);
R = blkdiag(R{:});

% Subtract mean from structure.
mu = 1 / N * S * ones(N, 1);
S = S - mu * ones(N, 1)';

% Project S on to basis.
S_sharp = k_reshape(S, 3);
[U, D, V] = svd(S_sharp, 'econ');
U = U(:, 1:K);
V = V(:, 1:K);
D = D(1:K, 1:K);
% Project to rank K.
S_sharp = U * D * V';
% Unreshape low-rank structure.
S_low_rank = k_unreshape(S_sharp, 3);
fprintf('Error to low rank = %g\n', ...
    norm(S - S_low_rank, 'fro') / norm(S, 'fro'));
S = S_low_rank;

% Find basis and coefficients such that S = kron(C, eye(3)) * B.
% Note these are not unique.
% B is 3N x K basis.
B = U;
% Unreshape to N x 3K.
B = k_unreshape(B, 3);
% C is F x K.
C = V * D;
P = R * kron(C, eye(3));

% Project.
W = R * S;

% Compute SVD of W.
[U, D, V] = svd(W, 'econ');
U = U(:, 1:3 * K);
V = V(:, 1:3 * K);
d = diag(D);
d = d(1:3 * K);

% Get initial factorization.
P_hat = 1 / sqrt(d(1)) * U * diag(sqrt(d));
B_hat = sqrt(d(1)) * diag(sqrt(d)) * V';

% Solve for corrective matrix.
Rs_hat = find_rotations(P_hat, 1e6);
%Rs_hat = find_rotations_dai(P_hat);
% Align cameras.
Rs_tilde = align_rotations(Rs_hat, Rs);
fprintf('Rotation error = %g\n', norm(Rs_tilde - Rs, 'fro') / norm(Rs, 'fro'));

% Convert to block-diagonal.
R_hat = mat2cell(Rs_hat, 2 * ones(F, 1), 3);
R_hat = cellfun(@sparse, R_hat, 'UniformOutput', false);
R_hat = blkdiag(R_hat{:});

fprintf('Any key to continue\n');
pause;

fprintf('Linear solution...\n');

[Q, C] = find_corrective_matrix(P_hat, R_hat);
S_linear = kron(C, eye(3)) * inv(Q) * B_hat;

fprintf('Reprojection error = %g\n', ...
    norm(W - R_hat * S_linear, 'fro') / norm(W, 'fro'));
fprintf('3D error = %g\n', ...
    norm(S - align_structure(S_linear, S), 'fro') / norm(S, 'fro'));

fprintf('Any key to continue\n');
pause;
fprintf('Nuclear norm solution...\n');

% Find structure.
S_nuclear = find_structure_affine_cameras(W, Rs_hat, true, ...
    struct(...
      'rho', 1, ...
      'mu', 10, ...
      'tau_incr', 2, ...
      'tau_decr', 2, ...
      'max_iter', 80, ...
      'epsilon_abs', 1e-3, ...
      'epsilon_rel', 1e-3, ...
      'min_rho_iter', 4));

fprintf('Reprojection error = %g\n', ...
    norm(W - R_hat * S_nuclear, 'fro') / norm(W, 'fro'));
fprintf('3D error = %g\n', ...
    norm(S - align_structure(S_nuclear, S), 'fro') / norm(S, 'fro'));

fprintf('Difference between solutions = %g\n', ...
    norm(S_linear - align_structure(S_nuclear, S_linear), 'fro') / ...
      sqrt(norm(S_linear, 'fro') * norm(S_linear, 'fro')));
