num_frames = 256;
num_sequences = 32;
downsample = 8;

omega_stddev = 5 * pi / 180;
scale_stddev = sqrt(2);
noise_stddevs = [0, 0.01, 1, 100];

num_noises = length(noise_stddevs);

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

K = 4;

% Project on to low rank.
for i = 1:num_sequences
  S = sequences(:, :, :, i);

  % Subtract 3D centroid from each frame.
  mu = mean(S, 2);
  S = bsxfun(@minus, S, mu);

  % [F, P, 3] -> [F, 3P]
  S = reshape(S, [num_frames, 3 * num_points]);
  % Project on to set of low rank matrices.
  S = project_rank(S, K);
  % [F, 3P] -> [F, P, 3]
  S = reshape(S, [num_frames, num_points, 3]);

  % Restore centroid.
  S = bsxfun(@plus, S, mu);

  sequences(:, :, :, i) = S;
end

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

% Solve reconstruction using each method.
residuals = zeros(num_sequences, num_noises, num_solvers);

for i = 1:num_sequences
  scene = scenes(i);

  for j = 1:num_noises
    noise_stddev = noise_stddevs(j);
    fprintf('Sequence %d\n', i);

    for k = 1:num_solvers
      solver = solvers(k);

      noise = noise_stddev * randn(num_frames, num_points, 2);
      solution = solver.solve(scene.projections + noise, K);

      residuals(i, j, k) = min_shape_error(scene.points, solution);
      residuals(i, j, k)
    end
  end
end

residuals = permute(residuals, [1, 3, 2]);

for i = 1:num_noises
  noise_stddev = noise_stddevs(i);
  figure;
  boxplot(residuals(:, :, i));
  hold on;
  plotSpread(residuals(:, :, i));
  hold off;
  grid on;
  axis([0, num_solvers + 1, 0, max(residuals(:)) * 1.05]);
  title(sprintf('Noise level %g', noise_stddev));
  set(gca, 'XTick', 1:num_solvers);
  set(gca, 'XTickLabel', {solvers.name});
  filename = sprintf('../figures/noise/%d-%g.eps', i, noise_stddev);
  print(filename, '-depsc2');
  unix(['epstopdf ', filename]);
  unix(['rm ', filename]);
end

residuals = permute(residuals, [1, 3, 2]);
