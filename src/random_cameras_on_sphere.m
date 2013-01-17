function cameras = random_cameras_on_sphere(n, radius, origin)
  cameras = zeros(3, 4, n);
  for t = 1:n
    cameras(:, :, t) = random_camera_on_sphere(radius, origin);
  end
end
