% Parameters:
% fig -- Figure handle
% points -- 3 x num_joints x num_frames matrix of joint positions.
% fps -- Upper limit on framerate to render at.
% lim -- Axis limits, can be empty.
% opts -- Additional arguments to plot().

function plot_movie(fig, points, fps, lim, opts)
  [d, N, F] = size(points);
  assert(d == 2 || d == 3, 'Data must be two- or three-dimensional');

  if isempty(lim)
    % [d, N, F] -> [NF, d]
    lim = axis_limits(reshape(points, [d, F * N])');
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
  animation.render = @(fig, i) render(fig, points, i, opts);
  animation.length = F;
  animation.fps = fps;

  animation.play();
end

function render(fig, points, i, opts)
  cla;
  plot_auto(gca(fig), points(:, :, i)', opts{:});
  title(sprintf('%d / %d', i, size(points, 3)));
end
