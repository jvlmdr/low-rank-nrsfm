function P = random_camera_on_sphere(radius, origin)
  % Random 3D unit vector.
  a = randn(3, 1);
  a = a / norm(a);

  % Another random 3D unit vector.
  b = randn(3, 1);
  b = b / norm(b);

  % Project b on to opposite hemisphere to a.

  % Measure component of b parallel to a.
  alpha = a' * b;
  % Project b perpendicular to a.
  b = b - alpha * a;
  % Subtract absolute component parallel to a.
  b = b - abs(alpha) * a;

  % Camera center.
  c = radius * a + origin;

  % Get direction from a to b. This will become camera 'look' vector.
  k = b - a;
  k = -k / norm(k);

  % Rotational freedom in choosing x and y.
  % Pick y to be most vertical.
  j = [0; 1; 0];
  i = cross(j, k);
  i = i / norm(i);
  j = cross(k, i);
  j = j / norm(j);
  R = [i, j, k]';

  % Build camera matrix.
  t = -R * c;
  P = [R, t];
end
