function points = load_mocap(input_file)
  % Load 3D points.
  data = load(input_file);
  points = data.S;

  num_frames = size(points, 1) / 3;
  num_points = size(points, 2);

  % [3F, N] -> [F, N, 3]
  points = shiftdim(reshape(points, [3, num_frames, num_points]), 1);

  % Swap back y and z.
  points = points(:, :, [1, 3, 2]);
end
