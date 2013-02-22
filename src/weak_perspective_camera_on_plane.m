% Parameters:
% theta -- scalar (radians, duh)
% scale -- scalar

function camera = weak_perspective_camera_on_plane(theta, scale)
  K = [scale 0 0 0; 0 scale 0 0; 0 0 0 1];
  R = roty(theta);
  t = zeros(3, 1);

  P = K * [R, t; zeros(1, 3), 1];

  camera = make_camera(P);
end
