num_sequences = numel(scenes);

num_joint = numel(solvers.find_both);
num_cameras = numel(solvers.find_cameras);
num_structure = numel(solvers.find_structure);
num_independent = num_cameras * num_structure;

errors_joint = nan(num_frames, num_joint, num_sequences);
errors_independent = nan(num_frames, num_independent, num_sequences);
benchmarks = nan(num_frames, num_sequences);

for i = 1:num_sequences
  S = scenes(i).points;

  for j = 1:num_joint
    X = solutions(i).joint(j).structure;
    for t = 1:num_frames
      errors_joint(t, j, i) = min_shape_error(S(:, :, t), X(:, :, t));
    end
  end

  for j = 1:num_independent
    X = solutions(i).independent(j).structure;
    if ~isempty(X)
      for t = 1:num_frames
        errors_independent(t, j, i) = min_shape_error(S(:, :, t), X(:, :, t));
      end
    end
  end

  projections = scenes(i).projections;
  X = [projections; zeros(1, num_points, num_frames)];
  for t = 1:num_frames
    benchmarks(t, i) = min_shape_error(S(:, :, t), X(:, :, t));
  end
end

errors_independent = reshape(errors_independent, num_frames, num_cameras, ...
    num_structure, num_sequences);

errors = struct(...
    'joint', errors_joint, ...
    'independent', errors_independent);

save('exact-rank-experiment-summary', 'errors', 'benchmarks');
