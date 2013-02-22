function joints = load_mocap(filename)
  % Load experiment data.
  mocap = load(filename);

  % Number of frames in experiment.
  num_frames = length(mocap.all_points);
  % Number of points in experiment.
  num_joints = length(mocap.segment_order);

  % Extract joints into one big matrix.
  joints = zeros(num_frames, num_joints, 3);
  for i = 1:num_frames
    joints(i, :, :) = mocap.all_points(i).pts';
  end
end
