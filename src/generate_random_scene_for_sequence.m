% Generates weak-perspective cameras with projections for a sequence.
%
% Parameters:
% points -- 3 x num_points x num_frames

function scene = generate_random_scene_for_sequence(points, omega_stddev, ...
    scale_stddev)
  num_frames = size(points, 3);
  num_points = size(points, 2);

  % Angular change in each frame.
  omegas = omega_stddev * randn(num_frames, 1);
  % Angle in each frame.
  thetas = cumsum(omegas) + rand() * 2 * pi;
  % Scale in each frame.
  scales = exp(log(scale_stddev) * randn(num_frames, 1));

  % Construct cameras.
  cameras = weak_perspective_cameras_on_plane(thetas, scales);

  % Project points.
  projections = project_points(cameras, points);

  scene = make_scene(points, cameras, projections);
end
