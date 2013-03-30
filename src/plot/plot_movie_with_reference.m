% plot_movie_with_reference(fig, points, fps, lim, opts)
%
% Parameters:
% fig -- Figure handle
% points -- 3 x num_joints x num_frames matrix of joint positions.
% ref_points -- 3 x num_joints x num_frames
% fps -- Upper limit on framerate to render at.
% lim -- Axis limits, can be empty.
% opts -- Additional arguments to plot().

function plot_movie_with_reference(fig, points, ref_points, fps, lim, opts, ...
    ref_opts)
  [d, N, F] = size(points);
  assert(d == 2 || d == 3, 'Data must be two- or three-dimensional');

  if isempty(lim)
    % [d, N, 2F] -> [2NF, d]
    lim = axis_limits(reshape(cat(3, points, ref_points), [d, 2 * F * N])');
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
  animation.render = @(fig, i) render(fig, points, ref_points, i, opts, ...
      ref_opts);
  animation.length = F;
  animation.fps = fps;

  animation.play();
end

function render(fig, points, ref_points, i, opts, ref_opts)
  F = size(points, 3);
  cla;
  plot_auto(gca(fig), points(:, :, i)', opts{:});
  hold on;
  plot_auto(gca(fig), ref_points(:, :, i)', ref_opts{:});
  title(sprintf('%d / %d', i, F));
end
