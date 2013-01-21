% Parameters:
% points -- num_frames x num_joints x 3 matrix of joint positions.
% ... -- Additional arguments to line().

function save_reconstruction_movie(fig, points, mask, colors, cameras, ...
    camera_size, lim, show_masked, plot_opts, masked_plot_opts, image_size, ...
    dpi, output_dir, movie_name)

  [num_frames, num_points, num_dims] = size(points);
  assert(num_dims == 3, 'Data must be three-dimensional');

  if isempty(lim)
    % Flatten out all points over all frames.
    range = reshape(points, [num_frames * num_points, num_dims]);

    % Get camera centers.
    centers = zeros(num_frames, 3);
    for t = 1:num_frames
      P = cameras(:, :, t);
      R = P(:, 1:3);
      d = P(:, 4);
      centers(t, :) = (-R' * d)';
    end

    lim = axis_limits([range; centers]);
  end

  init3(fig, lim);

  frame_dir = [output_dir, '/', movie_name];
  unix(['rm -rf ', frame_dir]);
  unix(['mkdir -p ', frame_dir]);

  frame_format = [frame_dir, '/%07d.png'];

  for t = 1:num_frames
    render_reconstruction_movie(fig, points, mask, colors, cameras, ...
        camera_size, t, show_masked, plot_opts, masked_plot_opts);
    drawnow;

    file = sprintf(frame_format, t);
    print_image(fig, image_size, dpi, file, {'-dpng'});
  end

  movie_file = [output_dir, '/', movie_name, '.mp4'];
  resize = sprintf('-vf scale=%u:%u', image_size(1), image_size(2));
  unix(['ffmpeg -sameq -y -i ', frame_format, ' ', resize, ' ', movie_file]);

  unix(['rm -rf ', frame_dir]);
end
