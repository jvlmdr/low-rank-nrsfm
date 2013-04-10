num_frames = 256;
num_sequences = 16;
downsample = 8;

omega_stddev = 5 * pi / 180;
scale_stddev = sqrt(2);

ranks = [2, 4, 8];

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

% Generate a camera for each sequence and project it.
scenes = generate_random_scene_for_all_sequences(sequences, omega_stddev, ...
    scale_stddev);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

find_cameras = [...
  make_solver(@find_rotations_rigid, 'Rigid', 'rigid'), ...
];

find_cameras_with_K = [...
  make_solver(@find_rotations_trace, 'Trace search', 'search'), ...
  make_solver(@find_rotations_dai, 'Dai et al', 'dai'), ...
];

find_structure = [...
  make_solver(...
    @(projections, rotations) find_structure_constrained_nuclear_norm(...
      projections, rotations, [], 2, 1e6, 1e-5, 1e-4, inf, true), ...
    'Nuclear', 'nuclear'), ...
];

find_structure_with_K = [...
  make_solver(@find_structure_nullspace, 'Nullspace', 'nullspace'), ...
];

solvers = struct(...
    'find_cameras', find_cameras, ...
    'find_cameras_with_K', find_cameras_with_K, ...
    'find_structure', find_structure, ...
    'find_structure_with_K', find_structure_with_K);

save('simple-experiment-setup', 'scenes', 'solvers', 'ranks');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trial = @(scene) { simple_experiment_trial(solvers, scene, ranks) };

if exist('pararrayfun', 'file')
  config = struct('h_cpu', '0:59:00', 'virtual_free', '1024M', ...
      'hostname', '!leffe*', 'matlab_optimisation', '1');
  solutions = pararrayfun(trial, scenes, numel(scenes), 'vp', config);
else
  warning('Could not find pararrayfun(), running in series.');
  solutions = arrayfun(trial, scenes(1));
end

% Repack.
num_solvers = numel(solutions{1});
solutions = reshape(cell2mat(solutions(:)), [size(solutions), num_solvers]);

save('simple-experiment-solutions', 'solutions');
