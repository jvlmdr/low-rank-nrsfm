% Generates weak-perspective cameras with projections for each sequence.
%
% Parameters:
% sequences -- num_frames x num_points x 3 x num_sequences

function scenes = generate_random_scene_for_all_sequences(sequences, ...
    omega_stddev, scale_stddev)
  num_frames = size(sequences, 1);
  num_points = size(sequences, 2);
  num_sequences = size(sequences, 4);

  % Convert to cell array of sequences.
  sequences = mat2cell(sequences, num_frames, num_points, 3, ...
      ones(num_sequences, 1));
  sequences = reshape(sequences, [num_sequences, 1]);

  scenes = cellfun(...
      @(sequence) generate_random_scene_for_sequence(sequence, omega_stddev, ...
        scale_stddev), ...
      sequences);
end
