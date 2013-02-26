clear;

num_frames = 256;
num_sequences = 32;
downsample = 8;

omega_stddev = 5 * pi / 180;
scale_stddev = sqrt(2);

ranks = [2, 4, 6, 8, 10];

solvers = [...
  make_solver(@(W, K) reconstruct_xiao_2004(W, K), 'Xiao (2004)', ...
    'xiao-2004'), ...
  make_solver(@(W, K) reconstruct_nuclear_norm(W, K), 'Nuclear norm', ...
    'nuclear-norm'), ...
  make_solver(@(W, K) reconstruct_nullspace(W, K), 'Nullspace', ...
    'nullspace'), ...
];

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

% Generate a camera for each sequence and project it.
scenes = generate_scene_for_each_sequence(sequences, omega_stddev, ...
    scale_stddev);

num_solvers = numel(solvers);
num_ranks = numel(ranks);

% Construct a K-specific solver for each method.
% Note this is a violation of the solver interface.
solvers_cross_ranks(num_solvers, num_ranks) = struct('solve', []);
for i = 1:num_solvers
  solver = solvers(i);
  for j = 1:num_ranks
    K = ranks(j);
    solvers_cross_ranks(i, j) = struct('solve', (@(W) solver.solve(W, K)));
  end
end

solve = @(projections, solver) struct('points', solver.solve(projections.W));
trial = @(projections) struct('solvers', ...
    arrayfun(@(solver) solve(projections, solver), solvers_cross_ranks));

% Extract projections in struct.
projections = struct('W', {scenes.projections});

if exist('pararrayfun', 'file')
  config = struct('virtual_free', '128M');
  results = pararrayfun(trial, projections, num_sequences, 'v', config);
else
  warning('Could not find pararrayfun(), running in series.');
  results = arrayfun(trial, projections);
end

% Solve reconstruction using each method.
residuals = zeros(num_sequences, num_solvers, num_ranks);

for i = 1:num_sequences
  scene = scenes(i);
  for j = 1:num_solvers
    solver = solvers(j);
    for k = 1:num_ranks
      K = ranks(k);
      solution = results(i).solvers(j, k);
      residuals(i, j, k) = min_shape_error(scene.points, solution.points);
    end
  end
end

save('../output/summary', 'residuals');
