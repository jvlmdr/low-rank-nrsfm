function render_reconstruction_movie(fig, points, mask, cameras, axis_length, t)
  cla;
  if isempty(mask)
    plot_auto(gca(fig), points(:, :, t), 'kx');
  else
    plot_auto(gca(fig), points(mask(t, :) ~= 0, :, t), 'kx');
    plot_auto(gca(fig), points(mask(t, :) == 0, :, t), 'yo');
  end
  plot_camera(gca(fig), cameras(:, :, t), axis_length);
  title(sprintf('%d / %d', t, size(points, 3)));
end
