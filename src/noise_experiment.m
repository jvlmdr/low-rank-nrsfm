num_frames = 200;
num_sequences = 20;
downsample = 8;

omega_stddev = 20 * pi / 180;
scale_stddev = 1.5;

% Load some mocap sequences.
if ~exist('../data/mocap-data.mat', 'file')
  sequences = random_mocap_sequences(NUM_FRAMES, NUM_SEQUENCES, DOWNSAMPLE, 42);
  save('../data/mocap-data', 'sequences');
else
  data = load('../data/mocap-data');
  sequences = data.sequences;
end

% In case these parameters differed in file which we loaded.
num_frames = size(sequences, 1);
num_points = size(sequences, 2);
num_sequences = size(sequences, 4);

% Generate a camera for each sequence and project it.
scenes = generate_scene_for_each_sequence(sequences, omega_stddev, ...
    scale_stddev);

solvers = [...
  make_solver(@(W, K) reconstruct_xiao_2004(W, K), 'Xiao (2004)', ...
    'xiao-2004'), ...
  make_solver(@(W, K) reconstruct_nuclear_norm(W, K), 'Nuclear norm', ...
    'nuclear-norm'), ...
  make_solver(@(W, K) reconstruct_nullspace(W, K), 'Nullspace', ...
    'nullspace'), ...
];

num_solvers = numel(solvers);

K = 8;
residuals = zeros(num_sequences, num_solvers);

% Solve reconstruction using each method.
for i = 1:num_sequences
  fprintf('Sequence %d\n', i);

  scene = scenes(i);

  for j = 1:num_solvers
    solver = solvers(j);

    solution = solver.solve(scene.projections, K);
    residuals(i, j) = min_shape_error(scene.points, solution);
  end
end
