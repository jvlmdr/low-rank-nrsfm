function render_camera_movie(fig, cameras, t, axis_length)
  ax = gca(fig);
  cla(ax);

  plot_camera(ax, cameras(:, :, t), axis_length);

  num_frames = size(cameras, 3);
  title(sprintf('%d / %d', t, num_frames));
end
