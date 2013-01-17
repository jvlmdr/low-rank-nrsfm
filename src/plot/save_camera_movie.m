% Parameters:
% cameras -- 3 x 4 x num_frames list of projection matrices

function save_camera_movie(fig, cameras, user_init, print_func, axis_length)

  num_frames = size(cameras, 3);

  % Get camera centers.
  centers = zeros(3, num_frames);
  for t = 1:num_frames
    P = cameras(:, :, t);
    R = P(:, 1:3);
    d = P(:, 4);
    centers(:, t) = -R' * d;
  end

  lim = axis_limits(centers');

  init3(fig, lim);

  if ~isempty(user_init)
    user_init(fig);
  end

  for t = 1:num_frames
    render_camera_movie(fig, cameras, t, axis_length);
    drawnow;
    print_func(fig, t);
  end
end
