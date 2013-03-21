clear;
close all;
rng('default');

K = 4;
scale_stddev = sqrt(2);
use_akhter_data = true;

if use_akhter_data
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
else
  % Load mocap sequence.
  data = load('../data/mocap-data.mat');
  F = size(data.sequences, 1);
  P = size(data.sequences, 2);
  points = data.sequences(:, :, :, 1);
  % Generate camera motion.
  scene = generate_scene_for_sequence(points, 20 * 180 / pi, scale_stddev);

  % Extract cameras.
  Rs = zeros(2, 3, F);
  for t = 1:F
    Rs(:, :, t) = scene.cameras(t).P(1:2, 1:3);
  end

  % [F, P, 3] -> [3, F, P]
  points = permute(points, [3, 1, 2]);
  % [3, F, P] -> [3F, P]
  S = points;
  S = reshape(points, [3 * F, P]);
end

% Build block-diagonal rotation matrix.
R = block_diagonal_cameras(Rs);

% Subtract centroid from structure.
mu = 1 / P * S * ones(P, 1);
S = S - mu * ones(P, 1)';

% Project S on to low-rank manifold.
S_sharp = k_reshape(S, 3);
S_sharp = project_rank(S_sharp, K);
S_low_rank = k_unreshape(S_sharp, 3);

%fprintf('Error to low rank = %g\n', ...
%    norm(S - S_low_rank, 'fro') / norm(S, 'fro'));
S = S_low_rank;

% Restore centroid.
S = S + mu * ones(P, 1)';

% [3F, P] -> [3, F, P] -> [3, P, F]
points = S;
points = reshape(points, [3, F, P]);
points = permute(points, [1, 3, 2]);

% Project.
W = R * S;

% Subtract mean from projections.
mu = 1 / P * W * ones(P, 1);
W = W - mu * ones(P, 1)';

% Build 2 x P x F matrix of projections.
% [2F, P] -> [2, F, P] -> [2, P, F]
projections = W;
projections = reshape(projections, [2, F, P]);
projections = permute(projections, [1, 3, 2]);

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

% [3F, P] -> [3, F, P] -> [3, P, F]
points_xiao = S_xiao;
points_xiao = reshape(points_xiao, [3, F, P]);
points_xiao = permute(points_xiao, [1, 3, 2]);
fprintf('3D error (Xiao 2004) = %g\n', min_shape_error(points, points_xiao));

%fprintf('Any key to continue\n');
%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified solution for cameras of Dai 2012

Rs_hat = find_rotations(M_hat, 1e6);
%Rs_hat = find_rotations_dai(M_hat);

% Convert to block-diagonal.
R_hat = block_diagonal_cameras(Rs_hat);

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

% [3F, P] -> [3, F, P] -> [3, P, F]
points_nuclear = S_nuclear;
points_nuclear = reshape(points_nuclear, [3, F, P]);
points_nuclear = permute(points_nuclear, [1, 3, 2]);;

fprintf('Reprojection error (nuclear) = %g\n', ...
    norm(W - R_hat * S_nuclear, 'fro') / norm(W, 'fro'));
fprintf('3D error (nuclear) = %g\n', min_shape_error(points, points_nuclear));

%fprintf('Any key to continue\n');
%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Our solution for structure using the nullspace

fprintf('Linear solution...\n');

[G, C] = find_corrective_matrix(M_hat, Rs_hat);
S_linear = kron(C, eye(3)) * inv(G) * B_hat;

% [3F, P] -> [3, F, P] -> [3, P, F]
points_linear = S_linear;
points_linear = reshape(points_linear, [3, F, P]);
points_linear = reshape(points_linear, [3, P, F]);

fprintf('Reprojection error (linear) = %g\n', ...
    norm(W - R_hat * S_linear, 'fro') / norm(W, 'fro'));
fprintf('3D error (linear) = %g\n', min_shape_error(points, points_linear));

%fprintf('Any key to continue\n');
%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for both structure and motion given estimate of R.

[Rs_nrsfm_nuclear, points_nrsfm_nuclear] = nrsfm_constrained_nuclear_norm(...
    projections, Rs_hat, 1, 1, 200, 10, 10, 10);

R_mat = block_diagonal_cameras(Rs_nrsfm_nuclear);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_mat = points_nrsfm_nuclear;
S_mat = permute(S_mat, [1, 3, 2]);
S_mat = reshape(S_mat, [3 * F, P]);

fprintf('Reprojection error (NRSFM nuclear) = %g\n', ...
    norm(W - R_mat * S_mat, 'fro') / norm(W, 'fro'));
fprintf('3D error (NRSFM nuclear) = %g\n', ...
    min_shape_error(points, points_nrsfm_nuclear));

%fprintf('Any key to continue\n');
%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for both structure and motion given estimate of R.

[Rs_nrsfm_rank, points_nrsfm_rank] = nrsfm_fixed_rank(projections, Rs_hat, K, 1, 1, ...
    200, 10, 10, 10);

R_mat = block_diagonal_cameras(Rs_nrsfm_rank);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_mat = points_nrsfm_rank;
S_mat = permute(S_mat, [1, 3, 2]);
S_mat = reshape(S_mat, [3 * F, P]);

fprintf('Reprojection error (NRSFM rank) = %g\n', ...
    norm(W - R_mat * S_mat, 'fro') / norm(W, 'fro'));
fprintf('3D error (NRSFM rank) = %g\n', ...
    min_shape_error(points, points_nrsfm_rank));

%fprintf('Any key to continue\n');
%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear refinement.

[Rs_nrsfm_nuclear_refined, points_nrsfm_nuclear_refined] = nrsfm_nonlinear(...
    projections, Rs_nrsfm_nuclear, points_nrsfm_nuclear, K, 1000, 1e-4);

R_mat = block_diagonal_cameras(Rs_nrsfm_nuclear_refined);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_mat = points_nrsfm_nuclear_refined;
S_mat = permute(S_mat, [1, 3, 2]);
S_mat = reshape(S_mat, [3 * F, P]);

fprintf('Reprojection error (non-linear) = %g\n', ...
    norm(W - R_mat * S_mat, 'fro') / norm(W, 'fro'));
fprintf('3D error (non-linear) = %g\n', ...
    min_shape_error(points, points_nrsfm_nuclear_refined));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear refinement.

[Rs_nrsfm_rank_refined, points_nrsfm_rank_refined] = nrsfm_nonlinear(...
    projections, Rs_nrsfm_rank, points_nrsfm_rank, K, 1000, 1e-4);

R_mat = block_diagonal_cameras(Rs_nrsfm_rank_refined);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_mat = points_nrsfm_rank_refined;
S_mat = permute(S_mat, [1, 3, 2]);
S_mat = reshape(S_mat, [3 * F, P]);

fprintf('Reprojection error (non-linear) = %g\n', ...
    norm(W - R_mat * S_mat, 'fro') / norm(W, 'fro'));
fprintf('3D error (non-linear) = %g\n', ...
    min_shape_error(points, points_nrsfm_rank_refined));
