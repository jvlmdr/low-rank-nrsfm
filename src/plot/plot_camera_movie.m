function plot_camera_movie(fig, C, lim, options)
  % Parameters:
  % fps -- Upper limit on framerate to render at.
  % ... -- Additional arguments to line().

  F = length(C);

  % Get camera centers.
  centers = zeros(F, 3);
  for t = 1:F
    P = C{t}.P;
    R = P(:, 1:3);
    d = P(:, 4);
    centers(t, :) = (-R' * d)';
  end

  if isempty(lim)
    lim = axis_limits(centers);
  end
  lim = lim([1, 2, 5, 6, 3, 4]);

  figure(fig);
  axis(lim);
  axis equal;
  axis manual;
  axis vis3d;
  hold on;
  set(gca(fig), 'YDir', 'reverse');
  grid on;

  animation = Animation(fig);
  animation.render = @(fig, t) render(fig, C, t);
  animation.length = F;
  animation.fps = fps;

  animation.play();
end

function render(fig, C, t)
  F = length(C);
  cla;
  plot_camera(gca(fig), C{t}.P);
  title(sprintf('%d / %d', t, F));
end
