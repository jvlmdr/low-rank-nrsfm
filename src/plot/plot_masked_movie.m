function plot_masked_movie(fig, points, mask, colors, lim, show_masked, ...
    plot_opts, masked_plot_opts, fps)
  % Parameters:
  % points -- num_frames x num_points x 3 matrix of joint positions.
  % mask -- num_frames x num_points visibility matrix
  % fps -- Upper limit on framerate to render at.
  % ... -- Additional arguments to line().

  [F, N, d] = size(points);
  assert(d == 2 || d == 3, 'Data must be two- or three-dimensional');

  if isempty(lim)
    lim = axis_limits(reshape(points, [F * N, d]));
  end
  if d == 3
    lim = lim([1, 2, 5, 6, 3, 4]);
  end

  figure(fig);
  axis(lim);
  axis equal;
  axis manual;
  if d == 3
    axis vis3d;
    set(gca(fig), 'YDir', 'reverse');
  end
  hold on;
  grid on;

  animation = Animation(fig);
  animation.render = @(fig, t) render_masked_movie(fig, points, mask, ...
      colors, t, show_masked, plot_opts, masked_plot_opts);
  animation.length = F;
  animation.fps = fps;

  animation.play();
end
