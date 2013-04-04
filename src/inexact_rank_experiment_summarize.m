num_frames = size(solutions(1).structure, 3);
num_points = size(solutions(1).structure, 2);
num_sequences = size(solutions, 1);
num_solvers = size(solutions, 3);

shape_errors = zeros(num_frames, num_sequences, num_solvers);
for i = 1:num_sequences
  S = scenes(i).points;

  for k = 1:num_solvers
    X = solutions(i, k).structure;

    for t = 1:num_frames
      shape_errors(t, i, k) = min_shape_error(S(:, :, t), X(:, :, t));
    end
  end
end

benchmark_errors = zeros(num_frames, num_sequences);
for i = 1:num_sequences
  S = scenes(i).points;

  projections = scenes(i).projections;
  X = [projections; zeros(1, num_points, num_frames)];

  for t = 1:num_frames
    benchmark_errors(t, i) = min_shape_error(S(:, :, t), X(:, :, t));
  end
end

save('inexact-experiment-summary', 'shape_errors', 'benchmark_errors');
