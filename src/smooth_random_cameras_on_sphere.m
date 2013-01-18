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
    % Pick x to be in direction of motion.
    if t == 1
      i = a(t + 1, :)' - a(t, :)';
    elseif t == n
      i = a(t, :)' - a(t - 1, :)';
    else
      i = a(t + 1, :)' - a(t - 1, :)';
    end
    j = cross(k, i);
    j = j / norm(j);
    i = cross(j, k);
    i = i / norm(i);
    R = [i, j, k]';

    % Build camera matrix.
    d = -R * c;
    P = [R, d];

    cameras(:, :, t) = P;
  end
end
