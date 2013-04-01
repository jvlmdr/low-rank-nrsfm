num_frames = size(solutions(1).structure, 3);
num_points = size(solutions(1).structure, 2);
num_sequences = size(solutions, 1);
num_noises = size(solutions, 2);
num_solvers = size(solutions, 3);

shape_errors = zeros(num_frames, num_sequences, num_noises, num_solvers);
for i = 1:num_sequences
  S = scenes(i).points;

  for j = 1:num_noises
    for k = 1:num_solvers
      X = solutions(i, j, k).structure;

      for t = 1:num_frames
        shape_errors(t, i, j, k) = min_shape_error(S(:, :, t), X(:, :, t));
      end
    end
  end
end

benchmark_errors = zeros(num_frames, num_sequences, num_noises);
for i = 1:num_sequences
  S = scenes(i).points;

  for j = 1:num_noises
    projections = noisy_scenes(i, j).projections;
    X = [projections; zeros(1, num_points, num_frames)];

    for t = 1:num_frames
      benchmark_errors(t, i, j) = min_shape_error(S(:, :, t), X(:, :, t));
    end
  end
end

save('noise-experiment-summary', 'shape_errors', 'benchmark_errors');
