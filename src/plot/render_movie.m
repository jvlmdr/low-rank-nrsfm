function render_movie(fig, points, colors, t, plot_opts)
  num_points = size(points, 2);

  % Shift the dimensions for easier access.
  points = shiftdim(points, 1);

  ax = gca(fig);
  cla(ax);

  for i = 1:num_points
    plot_auto(ax, points(i, :, t), plot_opts{:}, 'Color', colors(i, :));
  end
  title(sprintf('%d / %d', t, size(points, 3)));
end
