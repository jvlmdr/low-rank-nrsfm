clear all;
demo_setup;

% Using "ground truth" cameras (assuming that congealing finds transform
% yielding structure with optimal compaction).

structure = true_structure;
cameras = true_cameras;

S = k_reshape(structure_to_matrix(structure), 3);
S = project_rank(S, K);
structure = structure_from_matrix(k_unreshape(S, 3));

S = structure_to_matrix(structure);
R = block_diagonal_cameras(cameras);
W = R * S;
projections = projections_from_matrix(W);
structure_nullspace = find_structure_nullspace(projections, cameras, K);

keyboard;

%% Project true structure in to rank K volume.
%[basis, coeff] = factorize_structure(structure, K);
%structure_low_rank = compose_structure(basis, coeff);
%
%% Methods which only solve for structure, not cameras.
%% Find planar reconstruction.
%structure_planar = find_planar_structure(projections, cameras);
%
%structure_low_rank_refined = find_structure_low_rank_nonlinear(projections, ...
%    cameras, basis, coeff, 1000, 1e-4);
%structure_low_rank_refined_approx = find_structure_approx_low_rank_nonlinear(...
%    projections, structure, cameras, basis, coeff, 1, 1000, 1e-4);
%
%% Find min nuclear norm solution.
%structure_nuclear = find_structure_nuclear_norm_regularized(projections, [], ...
%    cameras, 1, 1e-6, 1.1, 1e6, 1e-4, 1e-3, inf, true);
%% Project nuclear norm solution in to rank K volume.
%[basis, coeff] = factorize_structure(structure_nuclear, K);
%structure_nuclear_low_rank = compose_structure(basis, coeff);
%structure_nuclear_refined = find_structure_low_rank_nonlinear(projections, ...
%    cameras, basis, coeff, 1000, 1e-4);
%structure_nuclear_refined_approx = find_structure_approx_low_rank_nonlinear(...
%    projections, structure_nuclear, cameras, basis, coeff, 1, 1000, 1e-4);
%
%fprintf('low_rank:                % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, structure_low_rank, cameras), ...
%    min_total_shape_error(unscaled_structure, structure_low_rank));
%fprintf('planar:                  % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, structure_planar, cameras), ...
%    min_total_shape_error(unscaled_structure, structure_planar));
%fprintf('low_rank_refined:        % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, structure_low_rank_refined, ...
%      cameras), ...
%    min_total_shape_error(unscaled_structure, structure_low_rank_refined));
%fprintf('low_rank_refined_approx: % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, ...
%      structure_low_rank_refined_approx, cameras), ...
%    min_total_shape_error(unscaled_structure, ...
%      structure_low_rank_refined_approx));
%fprintf('nuclear:                 % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, structure_nuclear, cameras), ...
%    min_total_shape_error(unscaled_structure, structure_nuclear));
%fprintf('nuclear_low_rank:        % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, structure_nuclear_low_rank, ...
%      cameras), ...
%    min_total_shape_error(unscaled_structure, structure_nuclear_low_rank));
%fprintf('nuclear_refined:         % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, structure_nuclear_refined, ...
%      cameras), ...
%    min_total_shape_error(unscaled_structure, structure_nuclear_refined));
%fprintf('nuclear_refined_approx:  % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, structure_nuclear_refined_approx, ...
%      cameras), ...
%    min_total_shape_error(unscaled_structure, ...
%      structure_nuclear_refined_approx));
%
%%% Refine nearest rank K matrix.
%%[structure_low_rank_refined_nuclear, cameras_low_rank_refined_nuclear] = ...
%%    nrsfm_nuclear_norm_regularized(projections, structure_low_rank, cameras, ...
%%    1, 1e-3, 1e-3, 1.1, 1e6, 1e-4, 1e-3, inf);
%%%[structure_low_rank_refined_nuclear, cameras_low_rank_refined_nuclear] = ...
%%%    nrsfm_nuclear_norm_regularized_alternation(projections, ...
%%%    structure_low_rank, cameras, 1, 1e-4, 100, [], 1.1, 1e6, 1e-4, 1e-3, inf);
%%[structure_low_rank_refined_nonlinear, cameras_low_rank_refined_nonlinear] = ...
%%    nrsfm_nonlinear(projections, cameras, basis, coeff, 1000, 1e-4);
%%% Refine nuclear norm solution.
%%[structure_nuclear_refined_nuclear, cameras_nuclear_refined_nuclear] = ...
%%    nrsfm_nuclear_norm_regularized(projections, structure_nuclear, cameras, ...
%%    1, 1e-3, 1e-3, 1.1, 1e6, 1e-4, 1e-3, inf);
%%%[structure_nuclear_refined_nuclear, cameras_nuclear_refined_nuclear] = ...
%%%    nrsfm_nuclear_norm_regularized_alternation(projections, ...
%%%    structure_low_rank, cameras, 1, 1e-4, 100, [], 1.1, 1e6, 1e-4, 1e-3, inf);
%%[basis, coeff] = factorize_structure(structure_nuclear, K);
%%structure_nuclear_low_rank = compose_structure(basis, coeff);
%%[structure_nuclear_refined_nonlinear, cameras_nuclear_refined_nonlinear] = ...
%%    nrsfm_nonlinear(projections, cameras, basis, coeff, 1000, 1e-4);
%%
%%fprintf('projection_error(planar) = %g\n', ...
%%    relative_projection_error(projections, structure_planar, cameras));
%%fprintf('projection_error(low_rank) = %g\n', ...
%%    relative_projection_error(projections, structure_low_rank, cameras));
%%fprintf('projection_error(low_rank_refined_nuclear) = %g\n', ...
%%    relative_projection_error(projections, ...
%%      structure_low_rank_refined_nuclear, cameras_low_rank_refined_nuclear));
%%fprintf('projection_error(low_rank_refined_nonlinear) = %g\n', ...
%%    relative_projection_error(projections, ...
%%      structure_low_rank_refined_nonlinear, ...
%%      cameras_low_rank_refined_nonlinear));
%%fprintf('projection_error(nuclear) = %g\n', ...
%%    relative_projection_error(projections, structure_nuclear, cameras));
%%fprintf('projection_error(nuclear_refined_nuclear) = %g\n', ...
%%    relative_projection_error(projections, ...
%%      structure_nuclear_refined_nuclear, cameras_nuclear_refined_nuclear));
%%fprintf('projection_error(nuclear_low_rank) = %g\n', ...
%%    relative_projection_error(projections, structure_nuclear_low_rank, ...
%%      cameras));
%%fprintf('projection_error(nuclear_refined_nonlinear) = %g\n', ...
%%    relative_projection_error(projections, ...
%%      structure_nuclear_refined_nonlinear, cameras_nuclear_refined_nonlinear));
%%
%%fprintf('shape_error(planar) = %g\n', ...
%%    min_total_shape_error(unscaled_structure, structure_planar));
%%fprintf('shape_error(low_rank) = %g\n', ...
%%    min_total_shape_error(unscaled_structure, structure_low_rank));
%%fprintf('shape_error(low_rank_refined_nuclear) = %g\n', ...
%%    min_total_shape_error(unscaled_structure, ...
%%      structure_low_rank_refined_nuclear));
%%fprintf('shape_error(low_rank_refined_nonlinear) = %g\n', ...
%%    min_total_shape_error(unscaled_structure, ...
%%      structure_low_rank_refined_nonlinear));
%%fprintf('shape_error(nuclear) = %g\n', ...
%%    min_total_shape_error(unscaled_structure, structure_nuclear));
%%fprintf('shape_error(nuclear_refined_nuclear) = %g\n', ...
%%    min_total_shape_error(unscaled_structure, ...
%%      structure_nuclear_refined_nuclear));
%%fprintf('shape_error(nuclear_low_rank) = %g\n', ...
%%    min_total_shape_error(unscaled_structure, structure_nuclear_low_rank));
%%fprintf('shape_error(nuclear_refined_nonlinear) = %g\n', ...
%%    min_total_shape_error(unscaled_structure, ...
%%      structure_nuclear_refined_nonlinear));
%
%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Solving rigid SfM\n');
rigid_cameras = find_rotations_rigid(projections);

fprintf('Solving cameras by Dai et al.''s method\n');
dai_cameras = find_rotations_dai(projections, K);

fprintf('Solving structure by Dai et al.''s method\n');
dai_structure = find_structure_dai(projections, dai_cameras, K);

% Estimate cameras.
fprintf('Solving cameras by our method\n');
cameras = find_rotations_trace(projections, K);

% Rotate "ground truth" structure into the reference frame of these cameras.
structure = zeros(3, P, F);
for t = 1:F
  U = cameras(:, :, t);
  U = [U; cross(U(1, :), U(2, :))];
  V = true_cameras(:, :, t);
  V = [V; cross(V(1, :), V(2, :))];
  % V S = U X = U (U' V S)
  structure(:, :, t) = U' * V * true_structure(:, :, t);
end

% Project true structure in to rank K volume.
[basis, coeff] = factorize_structure(structure, K);
structure_low_rank = compose_structure(basis, coeff);

% Methods which only solve for structure, not cameras.
% Find planar reconstruction.
structure_planar = find_planar_structure(projections, cameras);

structure_low_rank_refined = find_structure_low_rank_nonlinear(projections, ...
    cameras, basis, coeff, 1000, 1e-4);
structure_low_rank_refined_approx = find_structure_approx_low_rank_nonlinear(...
    projections, structure, cameras, basis, coeff, 1, 1000, 1e-4);

% Find min nuclear norm solution.
structure_nuclear = find_structure_constrained_nuclear_norm(projections, ...
    cameras, 1e-6, 2, 1e6, 1e-4, 1e-3, inf, true);
% Project nuclear norm solution in to rank K volume.
[basis, coeff] = factorize_structure(structure_nuclear, K);
structure_nuclear_low_rank = compose_structure(basis, coeff);
structure_nuclear_refined = find_structure_low_rank_nonlinear(projections, ...
    cameras, basis, coeff, 1000, 1e-4);
lambda = 1e-3;
structure_nuclear_refined_approx = find_structure_approx_low_rank_nonlinear(...
    projections, structure_nuclear, cameras, basis, coeff, lambda, 1000, 1e-4);

fprintf('low_rank:                % 12.3g % 12.3g\n', ...
    relative_projection_error(projections, structure_low_rank, cameras), ...
    min_total_shape_error(unscaled_structure, structure_low_rank));
fprintf('planar:                  % 12.3g % 12.3g\n', ...
    relative_projection_error(projections, structure_planar, cameras), ...
    min_total_shape_error(unscaled_structure, structure_planar));
fprintf('low_rank_refined:        % 12.3g % 12.3g\n', ...
    relative_projection_error(projections, structure_low_rank_refined, ...
      cameras), ...
    min_total_shape_error(unscaled_structure, structure_low_rank_refined));
fprintf('low_rank_refined_approx: % 12.3g % 12.3g\n', ...
    relative_projection_error(projections, ...
      structure_low_rank_refined_approx, cameras), ...
    min_total_shape_error(unscaled_structure, ...
      structure_low_rank_refined_approx));
fprintf('nuclear:                 % 12.3g % 12.3g\n', ...
    relative_projection_error(projections, structure_nuclear, cameras), ...
    min_total_shape_error(unscaled_structure, structure_nuclear));
fprintf('nuclear_low_rank:        % 12.3g % 12.3g\n', ...
    relative_projection_error(projections, structure_nuclear_low_rank, ...
      cameras), ...
    min_total_shape_error(unscaled_structure, structure_nuclear_low_rank));
fprintf('nuclear_refined:         % 12.3g % 12.3g\n', ...
    relative_projection_error(projections, structure_nuclear_refined, ...
      cameras), ...
    min_total_shape_error(unscaled_structure, structure_nuclear_refined));
fprintf('nuclear_refined_approx:  % 12.3g % 12.3g\n', ...
    relative_projection_error(projections, structure_nuclear_refined_approx, ...
      cameras), ...
    min_total_shape_error(unscaled_structure, ...
      structure_nuclear_refined_approx));

%% Project true structure in to rank K volume.
%[basis, coeff] = factorize_structure(structure, K);
%structure_low_rank = compose_structure(basis, coeff);
%
%% Methods which only solve for structure, not cameras.
%% Find planar reconstruction.
%structure_planar = find_planar_structure(projections, cameras);
%
%structure_low_rank_refined = find_structure_low_rank_nonlinear(projections, ...
%    cameras, basis, coeff, 1000, 1e-4);
%structure_low_rank_refined_approx = find_structure_approx_low_rank_nonlinear(...
%    projections, structure, cameras, basis, coeff, 1, 1000, 1e-4);
%
%% Find min nuclear norm solution.
%structure_nuclear = find_structure_nuclear_norm_regularized(projections, [], ...
%    cameras, 1, 1e-6, 1.1, 1e6, 1e-4, 1e-3, inf, true);
%% Project nuclear norm solution in to rank K volume.
%[basis, coeff] = factorize_structure(structure_nuclear, K);
%structure_nuclear_low_rank = compose_structure(basis, coeff);
%structure_nuclear_refined = find_structure_low_rank_nonlinear(projections, ...
%    cameras, basis, coeff, 1000, 1e-4);
%structure_nuclear_refined_approx = find_structure_approx_low_rank_nonlinear(...
%    projections, structure_nuclear, cameras, basis, coeff, 1, 1000, 1e-4);
%
%fprintf('low_rank:                % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, structure_low_rank, cameras), ...
%    min_total_shape_error(unscaled_structure, structure_low_rank));
%fprintf('planar:                  % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, structure_planar, cameras), ...
%    min_total_shape_error(unscaled_structure, structure_planar));
%fprintf('low_rank_refined:        % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, structure_low_rank_refined, ...
%      cameras), ...
%    min_total_shape_error(unscaled_structure, structure_low_rank_refined));
%fprintf('low_rank_refined_approx: % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, ...
%      structure_low_rank_refined_approx, cameras), ...
%    min_total_shape_error(unscaled_structure, ...
%      structure_low_rank_refined_approx));
%fprintf('nuclear:                 % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, structure_nuclear, cameras), ...
%    min_total_shape_error(unscaled_structure, structure_nuclear));
%fprintf('nuclear_low_rank:        % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, structure_nuclear_low_rank, ...
%      cameras), ...
%    min_total_shape_error(unscaled_structure, structure_nuclear_low_rank));
%fprintf('nuclear_refined:         % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, structure_nuclear_refined, ...
%      cameras), ...
%    min_total_shape_error(unscaled_structure, structure_nuclear_refined));
%fprintf('nuclear_refined_approx:  % 12.3g % 12.3g\n', ...
%    relative_projection_error(projections, structure_nuclear_refined_approx, ...
%      cameras), ...
%    min_total_shape_error(unscaled_structure, ...
%      structure_nuclear_refined_approx));

%% Find planar reconstruction.
%structure_planar = find_planar_structure(projections, cameras);
%% Project true structure on to rank K volume.
%[basis, coeff] = factorize_structure(structure, K);
%structure_low_rank = compose_structure(basis, coeff);
%% Refine nearest rank K matrix.
%[structure_low_rank_refined_nuclear, cameras_low_rank_refined_nuclear] = ...
%    nrsfm_nuclear_norm_regularized(projections, structure_nuclear, cameras, ...
%    1, 1e-3, 1e-3, 1.1, 1e6, 1e-4, 1e-3, inf);
%%[structure_low_rank_refined_nuclear, cameras_low_rank_refined_nuclear] = ...
%%    nrsfm_nuclear_norm_regularized_alternation(projections, ...
%%    structure_low_rank, cameras, 1, 1e-4, 100, [], 1.1, 1e6, 1e-4, 1e-3, inf);
%[structure_low_rank_refined_nonlinear, cameras_low_rank_refined_nonlinear] = ...
%    nrsfm_nonlinear(projections, cameras, basis, coeff, 1000, 1e-4);
%% Find min nuclear norm solution.
%structure_nuclear = find_structure_nuclear_norm_regularized(projections, [], ...
%    cameras, 1, 1e-6, 1.1, 1e6, 1e-4, 1e-3, inf, true);
%% Refine nuclear norm solution.
%[structure_nuclear_refined_nuclear, cameras_nuclear_refined_nuclear] = ...
%    nrsfm_nuclear_norm_regularized(projections, structure_nuclear, cameras, ...
%    1, 1e-3, 1e-3, 1.1, 1e6, 1e-4, 1e-3, inf);
%[basis, coeff] = factorize_structure(structure_nuclear, K);
%structure_nuclear_low_rank = compose_structure(basis, coeff);
%[structure_nuclear_refined_nonlinear, cameras_nuclear_refined_nonlinear] = ...
%    nrsfm_nonlinear(projections, cameras, basis, coeff, 1000, 1e-4);
%
%fprintf('projection_error(planar) = %g\n', ...
%    relative_projection_error(projections, structure_planar, cameras));
%fprintf('projection_error(low_rank) = %g\n', ...
%    relative_projection_error(projections, structure_low_rank, cameras));
%fprintf('projection_error(low_rank_refined_nuclear) = %g\n', ...
%    relative_projection_error(projections, ...
%      structure_low_rank_refined_nuclear, cameras_low_rank_refined_nuclear));
%fprintf('projection_error(low_rank_refined_nonlinear) = %g\n', ...
%    relative_projection_error(projections, ...
%      structure_low_rank_refined_nonlinear, ...
%      cameras_low_rank_refined_nonlinear));
%fprintf('projection_error(nuclear) = %g\n', ...
%    relative_projection_error(projections, structure_nuclear, cameras));
%fprintf('projection_error(nuclear_refined_nuclear) = %g\n', ...
%    relative_projection_error(projections, ...
%      structure_nuclear_refined_nuclear, cameras_nuclear_refined_nuclear));
%fprintf('projection_error(nuclear_low_rank) = %g\n', ...
%    relative_projection_error(projections, structure_nuclear_low_rank, ...
%      cameras));
%fprintf('projection_error(nuclear_refined_nonlinear) = %g\n', ...
%    relative_projection_error(projections, ...
%      structure_nuclear_refined_nonlinear, cameras_nuclear_refined_nonlinear));
%
%fprintf('shape_error(planar) = %g\n', ...
%    min_total_shape_error(unscaled_structure, structure_planar));
%fprintf('shape_error(low_rank) = %g\n', ...
%    min_total_shape_error(unscaled_structure, structure_low_rank));
%fprintf('shape_error(low_rank_refined_nuclear) = %g\n', ...
%    min_total_shape_error(unscaled_structure, ...
%      structure_low_rank_refined_nuclear));
%fprintf('shape_error(low_rank_refined_nonlinear) = %g\n', ...
%    min_total_shape_error(unscaled_structure, ...
%      structure_low_rank_refined_nonlinear));
%fprintf('shape_error(nuclear) = %g\n', ...
%    min_total_shape_error(unscaled_structure, structure_nuclear));
%fprintf('shape_error(nuclear_refined_nuclear) = %g\n', ...
%    min_total_shape_error(unscaled_structure, ...
%      structure_nuclear_refined_nuclear));
%fprintf('shape_error(nuclear_low_rank) = %g\n', ...
%    min_total_shape_error(unscaled_structure, structure_nuclear_low_rank));
%fprintf('shape_error(nuclear_refined_nonlinear) = %g\n', ...
%    min_total_shape_error(unscaled_structure, ...
%      structure_nuclear_refined_nonlinear));

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified solution for cameras of Dai 2012

estimated_cameras = find_rotations_trace(projections, K, 1e6);
%estimated_cameras = find_rotations_dai(M_hat);

% Get ground truth structure in estimated reference frame.
true_structure = zeros(3, P, F);
for t = 1:F
  U = estimated_rotations(:, :, t);
  U = [U; cross(U(1, :), U(2, :))];
  V = rotations(:, :, t);
  V = [V; cross(V(1, :), V(2, :))];
  % V S = U X = U (U' V S)
  true_structure(:, :, t) = U' * V * structure(:, :, t);
end
S_rel = structure_to_matrix(true_structure);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dai 2012 solution for structure

fprintf('Constrained nuclear norm solution...\n');

structure_hat = find_structure_constrained_nuclear_norm(projections, ...
    estimated_cameras, [], 400, 1e-6, 1.1, 1e6);
R_hat = block_diagonal_cameras(estimated_cameras);
S_hat = structure_to_matrix(structure_hat);

structure_planar = structure_from_matrix(pinv(full(R_hat)) * W);

fprintf('nuclear_norm(ground_truth) = %g\n', ...
    nuclear_norm(k_reshape(structure_to_matrix(true_structure), 3)));
fprintf('projection_error(ground_truth) = %g\n', ...
    norm(W - R_hat * S_rel, 'fro') / norm(W, 'fro'));

S_planar = structure_to_matrix(structure_planar);

fprintf('nuclear_norm(planar) = %g\n', ...
    nuclear_norm(k_reshape(structure_to_matrix(structure_planar), 3)));
fprintf('projection_error(planar) = %g\n', ...
    norm(W - R_hat * S_planar, 'fro') / norm(W, 'fro'));
fprintf('shape_error(planar) = %g\n', ...
    min_total_shape_error(unscaled_structure, structure_planar));


fprintf('nuclear_norm(solution) = %g\n', ...
    nuclear_norm(k_reshape(structure_to_matrix(structure_hat), 3)));
fprintf('projection_error(solution) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('shape_error(solution) = %g\n', ...
    min_total_shape_error(unscaled_structure, structure_hat));

structure_nuclear = structure_hat;
[basis_nuclear, coeff_nuclear] = factorize_structure(structure_nuclear, K);

keyboard;
%fprintf('Any key to continue\n');
%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Xiao 2004
%
%[structure_hat, rotations_hat] = nrsfm_basis_constraints(projections, K);
%
%R_hat = block_diagonal_cameras(rotations_hat);
%S_hat = structure_to_matrix(structure_hat);
%
%fprintf('Reprojection error (Xiao 2004) = %g\n', ...
%    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
%fprintf('3D error (Xiao 2004) = %g\n', ...
%    min_total_shape_error(unscaled_structure, structure_hat));
%
%%fprintf('Any key to continue\n');
%%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dai 2012 solution for structure, but with regularization not constraint.
%
%fprintf('Regularized nuclear norm solution...\n');
%
%lambda = 1e3;
%structure_hat = find_structure_nuclear_norm_regularized(projections, ...
%    estimated_cameras, lambda, [], 400, 1e-6, 1.1, 1e6);
%
%R_hat = block_diagonal_cameras(estimated_cameras);
%S_hat = structure_to_matrix(structure_hat);
%
%fprintf('projection_error(solution) = %g\n', ...
%    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
%fprintf('shape_error(solution) = %g\n', ...
%    min_total_shape_error(unscaled_structure, structure_hat));
%
%structure_nuclear_reg = structure_hat;
%
%keyboard;
%%fprintf('Any key to continue\n');
%%pause;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rank constraint by iterating Dai 2012 solution for structure
%
%fprintf('Rank problem by sweeping lambda...\n');
%
%[structure_hat, basis_hat, coeff_hat] = find_structure_nuclear_norm_sweep(...
%    projections, estimated_cameras, K, [], 200, 1e-6, 1.1, 1e6);
%
%R_hat = block_diagonal_cameras(estimated_cameras);
%S_hat = structure_to_matrix(structure_hat);
%
%fprintf('projection_error(solution) = %g\n', ...
%    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
%fprintf('shape_error(solution) = %g\n', ...
%    min_total_shape_error(unscaled_structure, structure_hat));
%
%keyboard;
%%fprintf('Any key to continue\n');
%%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear refinement.

for k = 1:5
  % Initialize with ground truth structure.
  R_hat = block_diagonal_cameras(estimated_cameras);
  [basis, coeff] = factorize_structure(true_structure, k);
  structure_low_rank = compose_structure(basis, coeff);
  S_low_rank = structure_to_matrix(structure_low_rank);

  [structure_hat, rotations_hat] = nrsfm_nonlinear(projections, ...
      estimated_cameras, basis, coeff, 200, 1e-4);

  fprintf('k = %d\n', k);
  fprintf('projection_error(ground_truth) = %g\n', ...
      norm(W - R_hat * S_rel, 'fro') / norm(W, 'fro'));

  fprintf('projection_error(low_rank) = %g\n', ...
      norm(W - R_hat * S_low_rank, 'fro') / norm(W, 'fro'));
  fprintf('shape_error(low_rank) = %g\n', ...
      min_total_shape_error(unscaled_structure, structure_low_rank));

  R_hat = block_diagonal_cameras(rotations_hat);
  S_hat = structure_to_matrix(structure_hat);
  fprintf('projection_error(non_linear) = %g\n', ...
      norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
  fprintf('shape_error(non_linear) = %g\n', ...
      min_total_shape_error(unscaled_structure, structure_hat));


  fprintf('projection_error(non_linear) = %g\n', ...
      norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
  fprintf('shape_error(solution) = %g\n', ...
      min_total_shape_error(unscaled_structure, structure_hat));

  x = align_all_shapes_similarity(unscaled_structure, structure_low_rank);
  y = align_all_shapes_similarity(unscaled_structure, structure_nuclear);
  fprintf('shape_diff(low_rank, nuclear_norm) = %g\n', ...
      norm(x(:) - y(:)) / sqrt(norm(x(:)) * norm(y(:))));

  x = align_all_shapes_similarity(unscaled_structure, structure_low_rank);
  y = align_all_shapes_similarity(unscaled_structure, structure_hat);
  fprintf('shape_diff(low_rank, non_linear) = %g\n', ...
      norm(x(:) - y(:)) / sqrt(norm(x(:)) * norm(y(:))));

  keyboard;
  %fprintf('Any key to continue\n');
  %pause;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear refinement.

[structure_hat, rotations_hat] = nrsfm_nonlinear(projections, ...
    rotations, basis_hat, coeff_hat, 1000, 1e-4);

R_hat = block_diagonal_cameras(rotations_hat);
S_hat = structure_to_matrix(structure_hat);

fprintf('Reprojection error (non-linear) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (non-linear) = %g\n', ...
    min_total_shape_error(unscaled_structure, structure_hat));

%fprintf('Any key to continue\n');
%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Our solution for structure using the nullspace

fprintf('Linear solution...\n');

structure_hat = find_structure_nullspace(projections, estimated_cameras, K);

R_hat = block_diagonal_cameras(estimated_cameras);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_hat = structure_hat;
S_hat = permute(S_hat, [1, 3, 2]);
S_hat = reshape(S_hat, [3 * F, P]);

fprintf('Reprojection error (linear) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (linear) = %g\n', ...
    min_total_shape_error(unscaled_structure, structure_hat));

%fprintf('Any key to continue\n');
%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve for both structure and motion given estimate of R.
%
%[Rs_nrsfm_nuclear, structure_nrsfm_nuclear] = nrsfm_constrained_nuclear_norm(...
%    projections, estimated_cameras, 1, 1, 200, 10, 10, 10);
%
%R_mat = block_diagonal_cameras(Rs_nrsfm_nuclear);
%% [3, P, F] -> [3, F, P] -> [3F, P]
%S_mat = structure_nrsfm_nuclear;
%S_mat = permute(S_mat, [1, 3, 2]);
%S_mat = reshape(S_mat, [3 * F, P]);
%
%fprintf('Reprojection error (NRSFM nuclear) = %g\n', ...
%    norm(W - R_mat * S_mat, 'fro') / norm(W, 'fro'));
%fprintf('3D error (NRSFM nuclear) = %g\n', ...
%    min_total_shape_error(unscaled_structure, structure_nrsfm_nuclear));
%
%%fprintf('Any key to continue\n');
%%pause;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nonlinear refinement of constrained nuclear norm ADMM solution for S and R.
%
%[Rs_refined_nrsfm_nuclear, structure_refined_nrsfm_nuclear] = ...
%    nrsfm_nonlinear(projections, Rs_nrsfm_nuclear, ...
%      structure_nrsfm_nuclear, K, 1000, 1e-4);
%
%R_mat = block_diagonal_cameras(Rs_refined_nrsfm_nuclear);
%% [3, P, F] -> [3, F, P] -> [3F, P]
%S_mat = structure_refined_nrsfm_nuclear;
%S_mat = permute(S_mat, [1, 3, 2]);
%S_mat = reshape(S_mat, [3 * F, P]);
%
%fprintf('Reprojection error (refined NRSFM nuclear) = %g\n', ...
%    norm(W - R_mat * S_mat, 'fro') / norm(W, 'fro'));
%fprintf('3D error (refined NRSFM nuclear) = %g\n', ...
%    min_total_shape_error(unscaled_structure, structure_refined_nrsfm_nuclear));
%
%%fprintf('Any key to continue\n');
%%pause;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve for both structure and motion given estimate of R.
%
%[Rs_nrsfm_rank, structure_nrsfm_rank] = nrsfm_fixed_rank(projections, ...
%    estimated_cameras, K, 1, 1, 200, 10, 10, 10);
%
%R_mat = block_diagonal_cameras(Rs_nrsfm_rank);
%% [3, P, F] -> [3, F, P] -> [3F, P]
%S_mat = structure_nrsfm_rank;
%S_mat = permute(S_mat, [1, 3, 2]);
%S_mat = reshape(S_mat, [3 * F, P]);
%
%fprintf('Reprojection error (NRSFM rank) = %g\n', ...
%    norm(W - R_mat * S_mat, 'fro') / norm(W, 'fro'));
%fprintf('3D error (NRSFM rank) = %g\n', ...
%    min_total_shape_error(unscaled_structure, structure_nrsfm_rank));
%
%%fprintf('Any key to continue\n');
%%pause;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nonlinear refinement of rank-constrained ADMM solution for S and R.
%
%[Rs_refined_nrsfm_rank, structure_refined_nrsfm_rank] = nrsfm_nonlinear(...
%    projections, Rs_nrsfm_rank, structure_nrsfm_rank, K, 1000, 1e-4);
%
%R_mat = block_diagonal_cameras(Rs_refined_nrsfm_rank);
%% [3, P, F] -> [3, F, P] -> [3F, P]
%S_mat = structure_refined_nrsfm_rank;
%S_mat = permute(S_mat, [1, 3, 2]);
%S_mat = reshape(S_mat, [3 * F, P]);
%
%fprintf('Reprojection error (refined NRSFM rank) = %g\n', ...
%    norm(W - R_mat * S_mat, 'fro') / norm(W, 'fro'));
%fprintf('3D error (refined NRSFM rank) = %g\n', ...
%    min_total_shape_error(unscaled_structure, structure_refined_nrsfm_rank));
%
%%fprintf('Any key to continue\n');
%%pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nullspace alternation, updating camera using motion matrix.
%
%clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;
%
%[structure_hat, rotations_hat] = nrsfm_nullspace_alternation_algebraic(...
%    projections, estimated_cameras, K, 40);
%
%R_hat = block_diagonal_cameras(rotations_hat);
%% [3, P, F] -> [3, F, P] -> [3F, P]
%S_hat = structure_hat;
%S_hat = permute(S_hat, [1, 3, 2]);
%S_hat = reshape(S_hat, [3 * F, P]);
%
%fprintf('Reprojection error (algebraic nullspace alternation) = %g\n', ...
%    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
%fprintf('3D error (algebraic nullspace alternation) = %g\n', ...
%    min_total_shape_error(unscaled_structure, structure_hat));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nullspace alternation, updating camera using structure.
%
%clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;
%
%[structure_hat, rotations_hat] = nrsfm_nullspace_alternation(...
%    projections, estimated_cameras, K, 40);
%
%R_hat = block_diagonal_cameras(rotations_hat);
%% [3, P, F] -> [3, F, P] -> [3F, P]
%S_hat = structure_hat;
%S_hat = permute(S_hat, [1, 3, 2]);
%S_hat = reshape(S_hat, [3 * F, P]);
%
%fprintf('Reprojection error (nullspace alternation) = %g\n', ...
%    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
%fprintf('3D error (nullspace alternation) = %g\n', ...
%    min_total_shape_error(unscaled_structure, structure_hat));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simple alternation, initialized using nullspace method.
%
%clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;
%
%[~, basis_hat] = find_structure_nullspace(projections, ...
%    estimated_cameras, K);
%[structure_hat, rotations_hat] = nrsfm_alternation(projections, ...
%    estimated_cameras, basis_hat, 80);
%
%R_hat = block_diagonal_cameras(rotations_hat);
%% [3, P, F] -> [3, F, P] -> [3F, P]
%S_hat = structure_hat;
%S_hat = permute(S_hat, [1, 3, 2]);
%S_hat = reshape(S_hat, [3 * F, P]);
%
%fprintf('Reprojection error (homogeneous alternation) = %g\n', ...
%    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
%fprintf('3D error (homogeneous alternation) = %g\n', ...
%    min_total_shape_error(unscaled_structure, structure_hat));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Minimize projection error regularized by nuclear norm.
%
%clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;
%
%lambda = 1;
%[structure_hat, rotations_hat] = nrsfm_nuclear_norm_regularized(...
%    projections, structure_nuclear_reg, estimated_cameras, lambda, 1, ...
%    200, 10, 10, 10);
%
%R_hat = block_diagonal_cameras(rotations_hat);
%% [3, P, F] -> [3, F, P] -> [3F, P]
%S_hat = structure_hat;
%S_hat = permute(S_hat, [1, 3, 2]);
%S_hat = reshape(S_hat, [3 * F, P]);
%
%fprintf('Reprojection error (nuclear norm regularizer) = %g\n', ...
%    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
%fprintf('3D error (nuclear norm regularizer) = %g\n', ...
%    min_total_shape_error(unscaled_structure, structure_hat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BALM

clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;

[structure_hat, rotations_hat] = nrsfm_balm_approximate(projections, ...
    estimated_cameras, coeff_nuclear, 1, 80, 10, 10, 10);

R_hat = block_diagonal_cameras(rotations_hat);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_hat = structure_hat;
S_hat = permute(S_hat, [1, 3, 2]);
S_hat = reshape(S_hat, [3 * F, P]);

fprintf('Reprojection error (BALM) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (BALM) = %g\n', ...
    min_total_shape_error(unscaled_structure, structure_hat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternation with metric projections

clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;

[structure_hat, rotations_hat] = nrsfm_metric_projections(projections, ...
    estimated_cameras, coeff_nuclear, 40);

R_hat = block_diagonal_cameras(rotations_hat);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_hat = structure_hat;
S_hat = permute(S_hat, [1, 3, 2]);
S_hat = reshape(S_hat, [3 * F, P]);

fprintf('Reprojection error (Alternation with metric projections) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (Alternation with metric projections) = %g\n', ...
    min_total_shape_error(unscaled_structure, structure_hat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BALM with metric projections

clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;

[structure_hat, rotations_hat] = nrsfm_balm_metric_projections(...
    projections, estimated_cameras, coeff_nuclear, 1, 40, 10, 10, 10);

R_hat = block_diagonal_cameras(rotations_hat);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_hat = structure_hat;
S_hat = permute(S_hat, [1, 3, 2]);
S_hat = reshape(S_hat, [3 * F, P]);

fprintf('Reprojection error (BALM with metric projections) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (BALM with metric projections) = %g\n', ...
    min_total_shape_error(unscaled_structure, structure_hat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternation using homogeneous problem

clear rotations_hat structure_hat basis_hat coeff_hat R_hat S_hat;

[structure_hat, rotations_hat] = nrsfm_homogeneous_alternation(...
    projections, estimated_cameras, basis_nuclear, 40);

R_hat = block_diagonal_cameras(rotations_hat);
% [3, P, F] -> [3, F, P] -> [3F, P]
S_hat = structure_hat;
S_hat = permute(S_hat, [1, 3, 2]);
S_hat = reshape(S_hat, [3 * F, P]);

fprintf('Reprojection error (homogeneous alternation) = %g\n', ...
    norm(W - R_hat * S_hat, 'fro') / norm(W, 'fro'));
fprintf('3D error (homogeneous alternation) = %g\n', ...
    min_total_shape_error(unscaled_structure, structure_hat));
