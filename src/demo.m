clear;
close all;
rng('default');

K = 4;
scale_stddev = sqrt(2);
use_akhter_data = false;

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

  % [F, P, 3] -> [3, F, P] -> [3F, P]
  S = points;
  S = permute(S, [3, 1, 2]);
  S = reshape(S, [3 * F, P]);
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

% Points with centroid removed for calculating error.
points_tilde = bsxfun(@minus, points, mean(points, 2));

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

[structure_hat, rotations_hat] = nrsfm_basis_constraints(projections, K);

R_hat = block_diagonal_cameras(rotations_hat);
S_hat = structure_to_matrix(structure_hat);

fprintf('Reprojection error (Xiao 2004) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (Xiao 2004) = %g\n', ...
    min_shape_error(points_tilde, structure_hat));

%fprintf('Any key to continue\n');
%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified solution for cameras of Dai 2012

rotations_trace = find_rotations(M_hat, 1e6);
%rotations_trace = find_rotations_dai(M_hat);

% Convert to block-diagonal.
R_hat = block_diagonal_cameras(rotations_trace);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dai 2012 solution for structure

fprintf('Nuclear norm solution...\n');

% Find structure.
S_nuclear = find_structure_affine_cameras(W, rotations_trace, true, ...
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
fprintf('3D error (nuclear) = %g\n', ...
    min_shape_error(points_tilde, points_nuclear));

%fprintf('Any key to continue\n');
%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear refinement of constrained nuclear norm solution for S.

%[Rs_refined_nuclear, points_refined_nuclear] = nrsfm_nonlinear(projections, ...
%    rotations_trace, points_nuclear, K, 1000, 1e-4);
%
%R_mat = block_diagonal_cameras(Rs_refined_nuclear);
%% [3, P, F] -> [3, F, P] -> [3F, P]
%S_mat = points_refined_nuclear;
%S_mat = permute(S_mat, [1, 3, 2]);
%S_mat = reshape(S_mat, [3 * F, P]);
%
%fprintf('Reprojection error (refined nuclear) = %g\n', ...
%    norm(W - R_mat * S_mat, 'fro') / norm(W, 'fro'));
%fprintf('3D error (refined nuclear) = %g\n', ...
%    min_shape_error(points_tilde, points_refined_nuclear));
%
%%fprintf('Any key to continue\n');
%%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Our solution for structure using the nullspace

fprintf('Linear solution...\n');

structure_hat = find_structure_nullspace(projections, rotations_trace, K);

R_hat = block_diagonal_cameras(rotations_trace);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_hat = structure_hat;
S_hat = permute(S_hat, [1, 3, 2]);
S_hat = reshape(S_hat, [3 * F, P]);

fprintf('Reprojection error (linear) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (linear) = %g\n', ...
    min_shape_error(points_tilde, structure_hat));

%fprintf('Any key to continue\n');
%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve for both structure and motion given estimate of R.
%
%[Rs_nrsfm_nuclear, points_nrsfm_nuclear] = nrsfm_constrained_nuclear_norm(...
%    projections, rotations_trace, 1, 1, 200, 10, 10, 10);
%
%R_mat = block_diagonal_cameras(Rs_nrsfm_nuclear);
%% [3, P, F] -> [3, F, P] -> [3F, P]
%S_mat = points_nrsfm_nuclear;
%S_mat = permute(S_mat, [1, 3, 2]);
%S_mat = reshape(S_mat, [3 * F, P]);
%
%fprintf('Reprojection error (NRSFM nuclear) = %g\n', ...
%    norm(W - R_mat * S_mat, 'fro') / norm(W, 'fro'));
%fprintf('3D error (NRSFM nuclear) = %g\n', ...
%    min_shape_error(points_tilde, points_nrsfm_nuclear));
%
%%fprintf('Any key to continue\n');
%%pause;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nonlinear refinement of constrained nuclear norm ADMM solution for S and R.
%
%[Rs_refined_nrsfm_nuclear, points_refined_nrsfm_nuclear] = nrsfm_nonlinear(...
%    projections, Rs_nrsfm_nuclear, points_nrsfm_nuclear, K, 1000, 1e-4);
%
%R_mat = block_diagonal_cameras(Rs_refined_nrsfm_nuclear);
%% [3, P, F] -> [3, F, P] -> [3F, P]
%S_mat = points_refined_nrsfm_nuclear;
%S_mat = permute(S_mat, [1, 3, 2]);
%S_mat = reshape(S_mat, [3 * F, P]);
%
%fprintf('Reprojection error (refined NRSFM nuclear) = %g\n', ...
%    norm(W - R_mat * S_mat, 'fro') / norm(W, 'fro'));
%fprintf('3D error (refined NRSFM nuclear) = %g\n', ...
%    min_shape_error(points_tilde, points_refined_nrsfm_nuclear));
%
%%fprintf('Any key to continue\n');
%%pause;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve for both structure and motion given estimate of R.
%
%[Rs_nrsfm_rank, points_nrsfm_rank] = nrsfm_fixed_rank(projections, ...
%    rotations_trace, K, 1, 1, 200, 10, 10, 10);
%
%R_mat = block_diagonal_cameras(Rs_nrsfm_rank);
%% [3, P, F] -> [3, F, P] -> [3F, P]
%S_mat = points_nrsfm_rank;
%S_mat = permute(S_mat, [1, 3, 2]);
%S_mat = reshape(S_mat, [3 * F, P]);
%
%fprintf('Reprojection error (NRSFM rank) = %g\n', ...
%    norm(W - R_mat * S_mat, 'fro') / norm(W, 'fro'));
%fprintf('3D error (NRSFM rank) = %g\n', ...
%    min_shape_error(points_tilde, points_nrsfm_rank));
%
%%fprintf('Any key to continue\n');
%%pause;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nonlinear refinement of rank-constrained ADMM solution for S and R.
%
%[Rs_refined_nrsfm_rank, points_refined_nrsfm_rank] = nrsfm_nonlinear(...
%    projections, Rs_nrsfm_rank, points_nrsfm_rank, K, 1000, 1e-4);
%
%R_mat = block_diagonal_cameras(Rs_refined_nrsfm_rank);
%% [3, P, F] -> [3, F, P] -> [3F, P]
%S_mat = points_refined_nrsfm_rank;
%S_mat = permute(S_mat, [1, 3, 2]);
%S_mat = reshape(S_mat, [3 * F, P]);
%
%fprintf('Reprojection error (refined NRSFM rank) = %g\n', ...
%    norm(W - R_mat * S_mat, 'fro') / norm(W, 'fro'));
%fprintf('3D error (refined NRSFM rank) = %g\n', ...
%    min_shape_error(points_tilde, points_refined_nrsfm_rank));
%
%%fprintf('Any key to continue\n');
%%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nullspace alternation, updating camera using motion matrix.

clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;

[structure_hat, rotations_hat] = nrsfm_nullspace_alternation_algebraic(...
    projections, rotations_trace, K, 40);

R_hat = block_diagonal_cameras(rotations_hat);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_hat = structure_hat;
S_hat = permute(S_hat, [1, 3, 2]);
S_hat = reshape(S_hat, [3 * F, P]);

fprintf('Reprojection error (algebraic nullspace alternation) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (algebraic nullspace alternation) = %g\n', ...
    min_shape_error(points_tilde, structure_hat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nullspace alternation, updating camera using structure.

clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;

[structure_hat, rotations_hat] = nrsfm_nullspace_alternation(projections, ...
    rotations_trace, K, 40);

R_hat = block_diagonal_cameras(rotations_hat);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_hat = structure_hat;
S_hat = permute(S_hat, [1, 3, 2]);
S_hat = reshape(S_hat, [3 * F, P]);

fprintf('Reprojection error (nullspace alternation) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (nullspace alternation) = %g\n', ...
    min_shape_error(points_tilde, structure_hat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple alternation, initialized using nullspace method.

clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;

[~, basis_hat] = find_structure_nullspace(projections, rotations_trace, K);
[structure_hat, rotations_hat] = nrsfm_alternation(projections, ...
    rotations_trace, basis_hat, 80);

R_hat = block_diagonal_cameras(rotations_hat);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_hat = structure_hat;
S_hat = permute(S_hat, [1, 3, 2]);
S_hat = reshape(S_hat, [3 * F, P]);

fprintf('Reprojection error (homogeneous alternation) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (homogeneous alternation) = %g\n', ...
    min_shape_error(points_tilde, structure_hat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimize projection error regularized by nuclear norm.

clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;

lambda = 1e3;
[structure_hat, rotations_hat] = nrsfm_nuclear_norm_regularizer(projections, ...
    rotations_trace, lambda, 1, 200, 10, 10, 10);

R_hat = block_diagonal_cameras(rotations_hat);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_hat = structure_hat;
S_hat = permute(S_hat, [1, 3, 2]);
S_hat = reshape(S_hat, [3 * F, P]);

fprintf('Reprojection error (nuclear norm regularizer) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (nuclear norm regularizer) = %g\n', ...
    min_shape_error(points_tilde, structure_hat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BALM, initialized using nullspace method.

clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;

[~, ~, coeff_hat] = find_structure_nullspace(projections, rotations_trace, K);
[structure_hat, rotations_hat] = nrsfm_balm_approximate(projections, ...
    rotations_trace, coeff_hat, 1, 80, 10, 10, 10);

R_hat = block_diagonal_cameras(rotations_hat);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_hat = structure_hat;
S_hat = permute(S_hat, [1, 3, 2]);
S_hat = reshape(S_hat, [3 * F, P]);

fprintf('Reprojection error (BALM) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (BALM) = %g\n', min_shape_error(points_tilde, structure_hat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternation with metric projections, initialized using nullspace method.

clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;

[~, ~, coeff_hat] = find_structure_nullspace(projections, rotations_trace, K);
[structure_hat, rotations_hat] = nrsfm_metric_projections(projections, ...
    rotations_trace, coeff_hat, 40);

R_hat = block_diagonal_cameras(rotations_hat);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_hat = structure_hat;
S_hat = permute(S_hat, [1, 3, 2]);
S_hat = reshape(S_hat, [3 * F, P]);

fprintf('Reprojection error (Alternation with metric projections) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (Alternation with metric projections) = %g\n', ...
    min_shape_error(points_tilde, structure_hat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BALM with metric projections, initialized using nullspace method.

clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;

[~, ~, coeff_hat] = find_structure_nullspace(projections, rotations_trace, K);
[structure_hat, rotations_hat] = nrsfm_balm_metric_projections(projections, ...
    rotations_trace, coeff_hat, 1, 40, 10, 10, 10);

R_hat = block_diagonal_cameras(rotations_hat);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_hat = structure_hat;
S_hat = permute(S_hat, [1, 3, 2]);
S_hat = reshape(S_hat, [3 * F, P]);

fprintf('Reprojection error (BALM with metric projections) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (BALM with metric projections) = %g\n', ...
    min_shape_error(points_tilde, structure_hat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternation using homogeneous problem, initialized using nullspace method.

clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;

[~, basis_hat] = find_structure_nullspace(projections, rotations_trace, K);
[structure_hat, rotations_hat] = nrsfm_homogeneous_alternation(projections, ...
    rotations_trace, basis_hat, 40);

R_hat = block_diagonal_cameras(rotations_hat);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_hat = structure_hat;
S_hat = permute(S_hat, [1, 3, 2]);
S_hat = reshape(S_hat, [3 * F, P]);

fprintf('Reprojection error (homogeneous alternation) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (homogeneous alternation) = %g\n', ...
    min_shape_error(points_tilde, structure_hat));
