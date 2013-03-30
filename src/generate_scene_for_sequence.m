% Generates weak-perspective cameras with projections for a sequence.
%
% Parameters:
% points -- 3 x num_points x num_frames

function scene = generate_scene_for_sequence(points, thetas, scales)
  % Construct cameras.
  cameras = weak_perspective_cameras_on_plane(thetas, scales);

  % Project points.
  projections = project_points(cameras, points);

  scene = make_scene(points, cameras, projections);
end
