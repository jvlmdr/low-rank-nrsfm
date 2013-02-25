clear;
close all;
rng('default');

K = 8;
scale_stddev = sqrt(2);

%% Load mocap sequence.
%data = load('../data/mocap-data.mat');
%F = size(data.sequences, 1);
%P = size(data.sequences, 2);
%points = data.sequences(:, :, :, 1);
%% Generate camera motion.
%scene = generate_scene_for_sequence(points, 20 * 180 / pi, 1);
%
%% [F, P, 3] -> [3, F, P] -> [3F, P]
%S = points;
%S = permute(S, [3, 1, 2]);
%S = reshape(S, [3 * F, P]);

% Load mocap sequence.
mocap_file = '../data/akhter-2008/yoga.mat';
load(mocap_file, 'S', 'Rs');
F = size(Rs, 1) / 2;
P = size(S, 2);

% [2F, 3] -> [2, F, 3] -> [2, 3, F]
Rs = reshape(Rs, [2, F, 3]);
Rs = permute(Rs, [1, 3, 2]);

% Scale each frame.
scales = exp(log(scale_stddev) * randn(F, 1));
Rs = bsxfun(@times, Rs, reshape(scales, [1, 1, F]));

% Build block-diagonal rotation matrix.
R = block_diagonal_cameras(Rs);

% Subtract centroid from structure.
mu = 1 / P * S * ones(P, 1);
S = S - mu * ones(P, 1)';

% Project S on to basis.
S_sharp = k_reshape(S, 3);
S_sharp = project_rank(S_sharp, K);
S_low_rank = k_unreshape(S_sharp, 3);

%fprintf('Error to low rank = %g\n', ...
%    norm(S - S_low_rank, 'fro') / norm(S, 'fro'));
S = S_low_rank;

% Restore centroid.
S = S + mu * ones(P, 1)';

% [3F, P] -> [3, F, P] -> [F, P, 3]
points = S;
points = reshape(points, [3, F, P]);
points = permute(points, [2, 3, 1]);

%% Find basis and coefficients such that S = kron(C, eye(3)) * B.
%% Note these are not unique.
%% B is 3N x K basis.
%B = U;
%% Unreshape to P x 3K.
%B = k_unreshape(B, 3);
%% C is F x K.
%C = V * D;
%M = R * kron(C, eye(3));

% Project.
W = R * S;

% Subtract mean from projections.
mu = 1 / P * W * ones(P, 1);
W = W - mu * ones(P, 1)';

% Compute SVD of W.
[U, D, V] = svd(W, 'econ');
U = U(:, 1:3 * K);
V = V(:, 1:3 * K);
d = diag(D);
d = d(1:3 * K);

% Get initial factorization.
M_hat = 1 / sqrt(d(1)) * U * diag(sqrt(d));
B_hat = sqrt(d(1)) * diag(sqrt(d)) * V';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xiao 2004

subset = randperm(F);
subset = subset(1:K);

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
R_hat = block_diagonal_cameras(Rs_hat);

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

% [3F, P] -> [3, F, P] -> [F, P, 3]
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

% [3F, P] -> [3, F, P] -> [F, P, 3]
points_linear = S_linear;
points_linear = reshape(points_linear, [3, F, P]);
points_linear = permute(points_linear, [2, 3, 1]);

fprintf('Reprojection error (linear) = %g\n', ...
    norm(W - R_hat * S_linear, 'fro') / norm(W, 'fro'));
fprintf('3D error (linear) = %g\n', min_shape_error(points, points_linear));
