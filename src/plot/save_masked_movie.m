% Parameters:
% points -- num_frames x num_joints x 3 matrix of joint positions.
% lim -- Manual axis limits. [] for automatic bounding box.
% plot_opts -- Additional arguments to plot().

function save_masked_movie(fig, points, mask, colors, lim, show_masked, ...
    plot_opts, masked_plot_opts, image_size, dpi, output_dir, movie_name)

  [num_frames, num_points, num_dims] = size(points);
  assert(num_dims == 2 || num_dims == 3, ...
      'Data must be two- or three-dimensional');

  if isempty(lim)
    range = reshape(points, [num_frames * num_points, num_dims]);
    if ~show_masked
      range = range(mask(:) ~= 0, :);
    end
    lim = axis_limits(range);
  end

  % Set axis limits.
  if num_dims == 2
    init2(fig, lim);
  else
    init3(fig, lim);
  end

  frame_dir = [output_dir, '/', movie_name];
  unix(['rm -rf ', frame_dir]);
  unix(['mkdir -p ', frame_dir]);

  frame_format = [frame_dir, '/%07d.png'];

  for t = 1:num_frames
    render_masked_movie(fig, points, mask, colors, t, show_masked, ...
        plot_opts, masked_plot_opts);
    drawnow;

    file = sprintf(frame_format, t);
    print_image(fig, image_size, dpi, file, {'-dpng'});
  end

  movie_file = [output_dir, '/', movie_name, '.mp4'];
  resize = sprintf('-vf scale=%u:%u', image_size(1), image_size(2));
  unix(['ffmpeg -sameq -y -i ', frame_format, ' ', resize, ' ', movie_file]);

  unix(['rm -rf ', frame_dir]);
end
