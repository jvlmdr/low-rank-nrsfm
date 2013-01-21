function render_masked_movie(fig, points, mask, colors, t, show_masked, ...
    plot_opts, masked_plot_opts)
  % Shift the dimensions for easier access.
  points = shiftdim(points, 1);
  num_points = size(points, 1);

  ax = gca(fig);
  cla(ax);

  for i = 1:num_points
    if isempty(mask) || mask(t, i)
      plot_auto(ax, points(i, :, t), plot_opts{:}, 'Color', colors(i, :));
    else
      if show_masked
        % Take color colors to white.
        %alpha = 0.3;
        %color = alpha * colors(i, :) + (1 - alpha);
        plot_auto(ax, points(i, :, t), masked_plot_opts{:}, ...
            'Color', colors(i, :));
      end
    end
  end

  title(sprintf('%d / %d', t, size(points, 3)));
end
