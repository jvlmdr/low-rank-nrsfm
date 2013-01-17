function print_image(fig, pixel_size, resolution, file, opts)
  inch_size = pixel_size / resolution;
  set(fig, 'PaperUnits','inches');
  set(fig, 'PaperSize', inch_size)
  set(fig, 'PaperPositionMode', 'manual');
  set(fig, 'PaperPosition', [0, 0, inch_size]);

  res_opt = sprintf('-r%d', 2 * resolution);
  print(fig, file, res_opt, opts{:});

  width = pixel_size(1);
  height = pixel_size(2);
  command = 'convert %s -resize %dx%d %s';
  unix(sprintf(command, file, width, height, file));
end
