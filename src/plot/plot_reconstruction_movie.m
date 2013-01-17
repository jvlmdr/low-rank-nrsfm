% Parameters:
% points -- num_frames x num_joints x 3 matrix of joint positions.
% fps -- Upper limit on framerate to render at.
% ... -- Additional arguments to line().

function plot_reconstruction_movie(fig, points, mask, cameras, fps, lim, ...
    axis_length)

  [F, N, dim] = size(points);
  assert(dim == 3, 'Data must be three-dimensional');

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

  % Shift the dimensions for easier access.
  points = shiftdim(points, 1);

  animation = Animation(fig);
  animation.render = @(fig, t) render_reconstruction_movie(fig, points, ...
      mask, cameras, axis_length, t);
  animation.length = F;
  animation.fps = fps;

  animation.play();
end
