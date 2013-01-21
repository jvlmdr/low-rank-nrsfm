% Parameters:
% points -- num_frames x num_joints x 3 matrix of joint positions.
% lim -- Manual axis limits. [] for automatic bounding box.
% plot_opts -- Additional arguments to plot().

function save_movie(fig, points, colors, lim, plot_opts, image_size, dpi, ...
    output_dir, movie_name)
  [num_frames, num_points, dim] = size(points);
  assert(dim == 2 || dim == 3, 'Data must be two- or three-dimensional');

  if isempty(lim)
    lim = axis_limits(reshape(points, [num_frames * num_points, dim]));
  end

  % Set axis limits.
  if dim == 2
    init2(fig, lim);
  else
    init3(fig, lim);
  end

  frame_dir = [output_dir, '/', movie_name];
  unix(['rm -rf ', frame_dir]);
  unix(['mkdir -p ', frame_dir]);

  frame_format = [frame_dir, '/%07d.png'];

  for t = 1:num_frames
    render_movie(fig, points, colors, t, plot_opts);
    drawnow;

    file = sprintf(frame_format, t);
    print_image(fig, image_size, dpi, file, {'-dpng'});
  end

  movie_file = [output_dir, '/', movie_name, '.mp4'];
  unix(['ffmpeg -sameq -y -i ', frame_format, ' ', movie_file]);

  unix(['rm -rf ', frame_dir]);
end
