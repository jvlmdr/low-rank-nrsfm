% Parameters:
% points -- num_frames x num_joints x 3 matrix of joint positions.
% ... -- Additional arguments to line().

function save_reconstruction_movie(fig, points, mask, cameras, lim, ...
    axis_length, user_init, print_func)

  [F, N, dim] = size(points);
  assert(dim == 3, 'Data must be three-dimensional');

  ax = gca(fig);

  % Get camera centers.
  centers = zeros(F, 3);
  for t = 1:F
    P = cameras(:, :, t);
    R = P(:, 1:3);
    d = P(:, 4);
    centers(t, :) = (-R' * d)';
  end

  if isempty(lim)
    lim = axis_limits([reshape(points, [F * N, dim]); centers]);
  end

  init3(fig, lim);

  if ~isempty(user_init)
    user_init(fig);
  end

  % Shift the dimensions for easier access.
  points = shiftdim(points, 1);

  for t = 1:F
    render_reconstruction_movie(fig, points, mask, cameras, axis_length, t);
    drawnow;
    print_func(fig, t);
  end
end
