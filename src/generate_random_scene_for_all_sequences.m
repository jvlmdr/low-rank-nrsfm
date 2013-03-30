% Generates weak-perspective cameras with projections for each sequence.
%
% Parameters:
% sequences -- 3 x num_points x num_frames x num_sequences

function scenes = generate_random_scene_for_all_sequences(sequences, ...
    omega_stddev, scale_stddev)
  num_points = size(sequences, 2);
  num_frames = size(sequences, 3);
  num_sequences = size(sequences, 4);

  % Convert to cell array of sequences.
  sequences = mat2cell(sequences, 3, num_points, num_frames, ...
      ones(num_sequences, 1));
  sequences = reshape(sequences, [num_sequences, 1]);

  scenes = cellfun(...
      @(sequence) generate_random_scene_for_sequence(sequence, omega_stddev, ...
        scale_stddev), ...
      sequences);
end
