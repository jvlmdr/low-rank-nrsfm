function cameras = smooth_random_cameras_on_sphere(n, r1, omega1, r2, ...
    omega2, origin)
  cameras = zeros(3, 4, n);

  a = generate_trajectory_on_sphere(n, r1, omega1);
  b = generate_trajectory_on_sphere(n, r2, omega2);

  for t = 1:n
    % Camera center.
    c = a(t, :)' + origin;
    % Direction from a to b. This will become camera 'look' vector.
    k = (b(t, :) - a(t, :))';
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
    d = -R * c;
    P = [R, d];

    cameras(:, :, t) = P;
  end
end
