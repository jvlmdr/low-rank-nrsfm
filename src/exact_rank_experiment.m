num_frames = 256;
num_sequences = 16;
downsample = 8;

omega_stddev = 5 * pi / 180;
scale_stddev = sqrt(2);

K = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load some mocap sequences.
if ~exist('../data/mocap-data.mat', 'file')
  sequences = random_mocap_sequences(num_frames, num_sequences, downsample, 42);
  save('../data/mocap-data', 'sequences');
else
  data = load('../data/mocap-data');
  sequences = data.sequences;
end

% In case these parameters differed in file which we loaded.
num_frames = size(sequences, 1);
num_points = size(sequences, 2);
num_sequences = min(num_sequences, size(sequences, 4));
sequences = sequences(:, :, :, 1:num_sequences);
% [F, P, 3, n] -> [3, P, F, n]
sequences = permute(sequences, [3, 2, 1, 4]);

% Project on to low rank.
for i = 1:num_sequences
  S = sequences(:, :, :, i);
  % Subtract 3D centroid from each frame.
  mu = mean(S, 2);
  S = bsxfun(@minus, S, mu);
  % Project on to rank-K manifold.
  S = k_reshape(structure_to_matrix(S), 3);
  S = project_rank(S, K);
  S = structure_from_matrix(k_unreshape(S, 3));
  % Restore centroid.
  S = bsxfun(@plus, S, mu);
  sequences(:, :, :, i) = S;
end

% Generate a camera for each sequence and project it.
scenes = generate_random_scene_for_all_sequences(sequences, omega_stddev, ...
    scale_stddev);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

find_both = [...
  make_solver(@(projections) nrsfm_basis_constraints(projections, K), ...
    'Xiao', 'xiao'), ...
];

find_cameras = [...
  make_solver(@(projections, true_cameras) true_cameras, 'True', 'true'), ...
  make_solver(...
    @(projections, true_cameras) find_rotations_rigid(projections), 'Rigid', ...
    'rigid'), ...
  make_solver(...
    @(projections, true_cameras) find_rotations_dai(projections, K), ...
    'Dai', 'dai'), ...
  make_solver(...
    @(projections, true_cameras) find_rotations_trace(projections, K), ...
    'Trace', 'trace'), ...
];

find_structure = [...
  make_solver(...
    @(projections, cameras) find_structure_constrained_nuclear_norm(...
      projections, cameras, [], 2, 1e6, 1e-5, 1e-4, inf, true), ...
    'Nuclear', 'nuclear'), ...
  make_solver(...
    @(projections, cameras) find_structure_nullspace(projections, cameras, ...
      K), ...
    'Nullspace', 'null'), ...
];

solvers = struct(...
    'find_both', find_both, ...
    'find_cameras', find_cameras, ...
    'find_structure', find_structure);

save('exact-rank-experiment-setup', 'scenes', 'solvers', 'K');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trial = @(scene) { exact_rank_experiment_trial(solvers, scene) };

if exist('pararrayfun', 'file')
  config = struct('h_cpu', '2:00:00', 'virtual_free', '1024M', ...
      'hostname', '!leffe*', 'matlab_optimisation', '1');
  solutions = pararrayfun(trial, scenes, min(numel(scenes), 8), 'vp', config);
else
  warning('Could not find pararrayfun(), running in series.');
  solutions = arrayfun(trial, scenes(1));
end

% Repack.
num_solvers = numel(solutions{1});
solutions = reshape(cell2mat(solutions(:)), [size(solutions), num_solvers]);

save('exact-rank-experiment-solutions', 'solutions');
