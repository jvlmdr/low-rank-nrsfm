% Plot the scene.
movie_name = 'scene';
fig = figure();
save_reconstruction_movie(fig, points, [], colors, extrinsics, 0.2, [], ...
    true, {'x', 'LineWidth', 2}, {}, image_size, dpi, movie_dir, movie_name);
close(fig);

for i = 1:num_solvers
  % Plot the reconstruction.
  structure = solutions(i).points;
  movie_name = ['solution-', solvers(i).id];
  fig = figure();
  save_reconstruction_movie(fig, structure, [], colors, extrinsics, 0.2, [], ...
      true, {'x', 'LineWidth', 2}, {'o', 'LineWidth', 1}, image_size, dpi, ...
      movie_dir, movie_name);
  close(fig);
end
