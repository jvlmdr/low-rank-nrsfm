num_sequences = numel(scenes);
num_ranks = numel(ranks);
num_cameras_sans = numel(solvers.find_cameras);
num_cameras_with = numel(solvers.find_cameras_with_K);
num_structure_sans = numel(solvers.find_structure);
num_structure_with = numel(solvers.find_structure_with_K);

num_cameras = num_cameras_sans + num_cameras_with;
num_structure = num_structure_sans + num_structure_with;

num_sans_sans = num_cameras_sans * num_structure_sans;
num_sans_with = num_ranks * num_cameras_sans * num_structure_with;
num_with_sans = num_ranks * num_cameras_with * num_structure_sans;
num_with_with = num_ranks * num_cameras_with * num_structure_with;
errors_sans_sans = zeros(num_frames, num_sans_sans, num_sequences);
errors_sans_with = zeros(num_frames, num_sans_with, num_sequences);
errors_with_sans = zeros(num_frames, num_with_sans, num_sequences);
errors_with_with = zeros(num_frames, num_with_with, num_sequences);

benchmarks = zeros(num_frames, num_sequences);

for i = 1:num_sequences
  S = scenes(i).points;

  for j = 1:num_sans_sans
    X = solutions(i).sans_sans(j).structure;
    for t = 1:num_frames
      errors_sans_sans(t, j, i) = min_shape_error(S(:, :, t), X(:, :, t));
    end
  end

  for j = 1:num_sans_with
    X = solutions(i).sans_with(j).structure;
    for t = 1:num_frames
      errors_sans_with(t, j, i) = min_shape_error(S(:, :, t), X(:, :, t));
    end
  end

  for j = 1:num_with_sans
    X = solutions(i).with_sans(j).structure;
    for t = 1:num_frames
      errors_with_sans(t, j, i) = min_shape_error(S(:, :, t), X(:, :, t));
    end
  end

  for j = 1:num_with_with
    X = solutions(i).with_with(j).structure;
    for t = 1:num_frames
      errors_with_with(t, j, i) = min_shape_error(S(:, :, t), X(:, :, t));
    end
  end

  projections = scenes(i).projections;
  X = [projections; zeros(1, num_points, num_frames)];
  for t = 1:num_frames
    benchmarks(t, i) = min_shape_error(S(:, :, t), X(:, :, t));
  end
end

errors_sans_sans = reshape(errors_sans_sans, num_frames, num_cameras_sans, ...
    num_structure_sans, num_sequences);
errors_sans_with = reshape(errors_sans_with, num_frames, num_ranks, ...
    num_cameras_sans, num_structure_with, num_sequences);
errors_with_sans = reshape(errors_with_sans, num_frames, num_ranks, ...
    num_cameras_with, num_structure_sans, num_sequences);
errors_with_with = reshape(errors_with_with, num_frames, num_ranks, ...
    num_cameras_with, num_structure_with, num_sequences);

errors = struct(...
    'sans_sans', errors_sans_sans, ...
    'sans_with', errors_sans_with, ...
    'with_sans', errors_with_sans, ...
    'with_with', errors_with_with);

save('simple-experiment-summary', 'errors', 'benchmarks');
