NUM_FRAMES = 200;
NUM_SEQUENCES = 20;
DOWNSAMPLE = 8;

% Load some mocap sequences.
if ~exist('../data/mocap-data.mat', 'file')
  sequences = random_mocap_sequences(NUM_FRAMES, NUM_SEQUENCES, DOWNSAMPLE, 42);
  save('../data/mocap-data', '-struct', 'data');
else
  data = load('../data/mocap-data');
end

% Generate a camera trajectory for each sequence.

% Project each sequence.

% Solve reconstruction using each method.

% Compute error.
