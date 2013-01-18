function x = generate_trajectory_on_sphere(n, r, sigma)

  % Start at a random point on the sphere.
  y = randn(3, 1);
  y = y / norm(y);

  % Generate a 3D velocity.
  dq = sigma * generate_velocity(n - 1, 3);
  % Divide steps by radius to obtain angles in radians.
  dq = dq / r;

  x = y';
  for i = 1:n - 1
    p = x(i, :)';
    p = rotx(dq(i, 1)) * p;
    p = roty(dq(i, 2)) * p;
    p = rotz(dq(i, 3)) * p;

    % Ensure still lies on the sphere (bit paranoid...)
    x(i + 1, :) = p' / norm(p);
  end

  x = r * x;
end
