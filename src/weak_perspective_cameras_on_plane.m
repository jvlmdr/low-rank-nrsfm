% Parameters:
% theta -- num_frames x 1 (radians, duh)
% scale -- num_frames x 1

function cameras = weak_perspective_cameras_on_plane(thetas, scales)
  num_frames = length(thetas);

  cameras(num_frames, 1) = make_camera();

  for t = 1:num_frames
    cameras(t) = weak_perspective_camera_on_plane(thetas(t), scales(t));
  end
end
