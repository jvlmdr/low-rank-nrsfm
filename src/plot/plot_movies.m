function plot_movies(fig, movies, fps, opts)
  % Parameters:
  % movies -- rows x cols struct array, each with
  %   points -- num_frames x num_joints x 3 matrix of joint positions
  %   title
  % fps -- Upper limit on framerate to render at.
  % ... -- Additional arguments to line().

  [rows, cols] = size(movies);
  [num_frames, num_points, three] = size(movies(1).points);

  % Select active figure and clear.
  figure(fig);
  clf(fig);

  for i = 1:rows
    for j = 1:cols
      points = movies(i, j).points;

      % Calculate bounds for axes.
      lim = axis_limits(reshape(points, [num_frames * num_points, 3]));

      % Set up axes.
      subplot(rows, cols, (j - 1) * rows + i);
      %title(movies.title);
      axis(lim);
      axis equal;
      axis manual;
      axis vis3d;
      hold on;
    end
  end

  animation = Animation(fig);
  animation.render = @(fig, t) render(fig, movies, t, opts);
  animation.length = num_frames;
  animation.fps = fps;

  animation.play();
end

function render(fig, movies, t, opts)
  [rows, cols] = size(movies);

  figure(fig);

  for i = 1:rows
    for j = 1:cols
      points = movies(i, j).points;
      subplot(rows, cols, (j - 1) * rows + i);
      cla;
      points = shiftdim(points, 1);
      plot23(points(:, :, t), opts{:});
    end
  end
end
