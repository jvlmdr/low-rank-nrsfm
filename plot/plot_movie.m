function plot_movie(fig, points, fps, lim, opts)
  % Parameters:
  % points -- num_frames x num_joints x 3 matrix of joint positions.
  % fps -- Upper limit on framerate to render at.
  % ... -- Additional arguments to line().

  [F, N, d] = size(points);
  assert(d == 2 || d == 3, 'Data must be two- or three-dimensional');

  if isempty(lim)
    lim = axis_limits(reshape(points, [F * N, d]));
  end

  figure(fig);
  axis(lim);
  axis equal;
  axis manual;
  if d == 3
    axis vis3d;
  end
  hold on;

  % Shift the dimensions for easier access.
  points = shiftdim(points, 1);

  animation = Animation(fig);
  animation.render = @(fig, i) render(fig, points, i, opts);
  animation.length = F;
  animation.fps = fps;

  animation.play();
end

function render(fig, points, i, opts)
  cla;
  plot_auto(gca(fig), points(:, :, i), opts{:});
end
