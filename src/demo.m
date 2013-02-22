clear;
close all;
rng('default');

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
M = R * kron(C, eye(3));

% Project.
W = R * S;

% Compute SVD of W.
[U, D, V] = svd(W, 'econ');
U = U(:, 1:3 * K);
V = V(:, 1:3 * K);
d = diag(D);
d = d(1:3 * K);

% Get initial factorization.
M_hat = 1 / sqrt(d(1)) * U * diag(sqrt(d));
B_hat = sqrt(d(1)) * diag(sqrt(d)) * V';

subset = randperm(F);
subset = subset(1:K);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xiao 2004

% Solve for corrective matrix.
[G, Rs_hat, C_hat] = find_corrective_transform_xiao_2004_linear(M_hat, subset);
% Recover structure.
S_xiao = kron(C_hat, eye(3)) * inv(G) * B_hat;

% [3F, P] -> [3, F, P] -> [F, P, 3]
points_xiao = S_xiao;
points_xiao = reshape(points_xiao, [3, F, P]);
points_xiao = permute(points_xiao, [2, 3, 1]);
fprintf('3D error (Xiao 2004) = %g\n', min_shape_error(points, points_xiao));

fprintf('Any key to continue\n');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified solution for cameras of Dai 2012

Rs_hat = find_rotations(M_hat, 1e6);
%Rs_hat = find_rotations_dai(M_hat);

% Convert to block-diagonal.
R_hat = mat2cell(Rs_hat, 2 * ones(F, 1), 3);
R_hat = cellfun(@sparse, R_hat, 'UniformOutput', false);
R_hat = blkdiag(R_hat{:});

fprintf('Any key to continue\n');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dai 2012 solution for structure

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

points_nuclear = S_nuclear;
points_nuclear = reshape(points_nuclear, [3, F, P]);
points_nuclear = permute(points_nuclear, [2, 3, 1]);

fprintf('Reprojection error (nuclear) = %g\n', ...
    norm(W - R_hat * S_nuclear, 'fro') / norm(W, 'fro'));
fprintf('3D error (nuclear) = %g\n', min_shape_error(points, points_nuclear));

fprintf('Any key to continue\n');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Our solution for structure using the nullspace

fprintf('Linear solution...\n');

[G, C] = find_corrective_matrix(M_hat, R_hat);
S_linear = kron(C, eye(3)) * inv(G) * B_hat;

points_linear = S_linear;
points_linear = reshape(points_linear, [3, F, P]);
points_linear = permute(points_linear, [2, 3, 1]);

fprintf('Reprojection error (linear) = %g\n', ...
    norm(W - R_hat * S_linear, 'fro') / norm(W, 'fro'));
fprintf('3D error (linear) = %g\n', min_shape_error(points, points_linear));
