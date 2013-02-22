function [joints, sequence_files, offsets] = random_mocap_sequences(...
    desired_length, desired_sequences, downsample, seed)

  % Directory to load mocap data from.
  MOCAP_DIR = '/staticdata/cmu-mocap/';
  INDEX_FILE = 'all-point-paths';

  % Number of frames we need before downsampling!
  n = desired_length * downsample;

  % Load master mocap file.
  mocap = load([MOCAP_DIR, INDEX_FILE]);
  num_sequences = length(mocap.all_paths);

  % Load each sequence and find its length.
  num_frames = nan(1, num_sequences);
  for i = 1:num_sequences
    if mod(i, 100) == 0
      fprintf('Loading sequence %d\n', i);
    end
    filename = [MOCAP_DIR, mocap.all_paths(i).path];
    sequence_data = load(filename);
    num_frames(i) = length(sequence_data.all_points);
  end

  % Filter sequences by length.
  candidates = find(num_frames >= n);

  % Generate way too many random numbers (to ensure unique).
  stream = RandStream('mrg32k3a', 'Seed', seed);
  index = stream.randperm(length(candidates));
  % And trim down to the number required.
  index = index(1:desired_sequences);

  % Get index of each sequence in the final set.
  sequences = candidates(index);
  offsets = zeros(size(sequences));
  % Get number of points from last sequence.
  num_points = size(sequence_data.all_points(1).pts, 2);
  % Initialize structure to put sequence data in.
  joints = zeros(desired_length, num_points, 3, desired_sequences);
  sequence_files = cell(1, desired_sequences);

  % Load each sequence.
  for i = 1:desired_sequences
    sequence = sequences(i);

    % Load corresponding mocap file.
    rel_path = mocap.all_paths(sequence).path;
    filename = fullfile(MOCAP_DIR, rel_path);
    sequence_data = load_mocap(filename);

    % If the sequence is too long, pick some random section.
    extra = num_frames(sequence) - n;
    offset = randint(1, 1, extra + 1);

    % Extract points from each frame.
    frames = downsample * ((1:desired_length) - 1) + 1 + offset;
    joints(:, :, :, i) = sequence_data(frames, :, :);
    offsets(i) = offset;
    sequence_files{i} = rel_path;
  end
end
