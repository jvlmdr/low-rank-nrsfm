function points = load_mocap(input_file)
  % Load 3D points.
  data = load(input_file);

  num_frames = length(data.all_points);
  num_points = size(data.all_points(1).pts, 2);

  points = zeros(3, num_points, num_frames);
  for t = 1:num_frames
    points(:, :, t) = data.all_points(t).pts;
  end

  points = permute(points, [3, 2, 1]);
end
