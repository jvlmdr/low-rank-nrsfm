function save_movies(fig, movies, format, plot_opts, print_opts)
  % Parameters:
  % movies -- rows x cols struct array, each with
  %   points -- num_frames x num_joints x 3 matrix of joint positions
  %   title
  % plot_opts -- Additional arguments to plot().
  % print_opts -- Additional arguments to print().

  [rows, cols] = size(movies);
  [num_frames, num_points, three] = size(movies(1).points);

  for i = 1:rows
    for j = 1:cols
      points = movies(i, j).points;

      % Create axes.
      ax(i, j) = subplot(rows, cols, (i - 1) * cols + j);
      set(ax(i, j), 'FontSize', 16);

      % Set up axes.
      if isfield(movies(i, j), 'title')
        title(ax(i, j), movies(i, j).title);
      end

      % Calculate bounds for axes.
      lim = axis_limits(reshape(points, [num_frames * num_points, 3]));
      axis(ax(i, j), lim);

      axis(ax(i, j), 'equal');
      axis(ax(i, j), 'manual');
      axis(ax(i, j), 'vis3d');
      axis(ax(i, j), 'off');
      hold(ax(i, j), 'on');
    end
  end

  for t = 1:num_frames
    render(ax, movies, t, plot_opts);
    %drawnow;
    file = sprintf(format, t);
    print(fig, file, print_opts{:});
  end
end

function render(ax, movies, t, opts)
  [rows, cols] = size(movies);

  if mod(t, 10) == 0
    n = size(movies(1, 1).points, 1);
    fprintf('Frame %d / %d\n', t, n);
  end

  for i = 1:rows
    for j = 1:cols
      points = shiftdim(movies(i, j).points, 1);
      cla(ax(i, j));
      plot_auto(ax(i, j), points(:, :, t), opts{:});
    end
  end
end
