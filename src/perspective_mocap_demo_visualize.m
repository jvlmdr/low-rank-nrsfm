% Visualize input sequence with different missing data modes.

for j = 1:length(kinds_of_missing_data)
  index = kinds_of_missing_data(j).index;
  mask = kinds_of_missing_data(j).mask;
  missing_data_id = kinds_of_missing_data(j).id;

  % Plot the 3D input with patterns of missing data.
  movie_name = ['missing-data-', missing_data_id];
  fig = figure();
  save_masked_movie(fig, points(:, index, :), mask, colors, [], ...
      false, {'x', 'LineWidth', 2}, {'o', 'LineWidth', 1}, image_size, dpi, ...
      movie_dir, movie_name);
  close(fig);
end

for i = 1:length(camera_paths)
  camera_path_id = camera_paths(i).id;
  projections = camera_paths(i).projections;
  extrinsics = camera_paths(i).extrinsics;

%  % Plot the projections.
%  movie_name = ['projections-', camera_path_id];
%  fig = figure();
%  save_masked_movie(fig, projections(:, :, :), ones(num_frames, num_points), ...
%      colors, [], false, {'x', 'LineWidth', 2}, {'o', 'LineWidth', 1}, ...
%      image_size, dpi, movie_dir, movie_name);
%  close(fig);

  % Plot the sequence with the camera motion.
  structure = points;
  movie_name = ['scene-', camera_path_id];
  fig = figure();
  save_reconstruction_movie(fig, structure(:, index, :), [], colors, ...
      extrinsics, 0.2, [], true, {'x', 'LineWidth', 2}, {}, image_size, dpi, ...
      movie_dir, movie_name);
  close(fig);

  for j = 1:length(kinds_of_missing_data)
    index = kinds_of_missing_data(j).index;
    mask = kinds_of_missing_data(j).mask;
    missing_data_id = kinds_of_missing_data(j).id;

    % Plot the reconstruction.
    structure = solutions(i, j).points;
    movie_name = ['solution-', camera_path_id, '-', missing_data_id];
    fig = figure();
    save_reconstruction_movie(fig, structure(:, index, :), mask, colors, ...
        extrinsics, 0.2, [], true, {'x', 'LineWidth', 2}, ...
        {'o', 'LineWidth', 1}, image_size, dpi, movie_dir, movie_name);
    close(fig);
  end
end
