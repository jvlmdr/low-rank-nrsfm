% ADMM settings.
settings = struct(...
    'rho', 1, ...
    'mu', 10, ...
    'tau_incr', 2, ...
    'tau_decr', 2, ...
    'max_iter', 100, ...
    'epsilon_abs', 1e-3, ...
    'epsilon_rel', 1e-3);

solutions = struct('points', {});

num_camera_paths = size(sequences, 1);
num_sequences = size(sequences, 2);

% For each set of observations, reconstruct.
for i = 1:num_camera_paths
  cameras = camera_paths(i).cameras;

  for j = 1:num_sequences
    observations = sequences(i, j);
    tracks = observations.tracks;
    num_tracks = length(tracks);
    mask = kinds_of_missing_data(j).mask;
    index = kinds_of_missing_data(j).index;

    % Visualize projections.
    if interactive
      % Plot projections.
      fig = figure();
      plot_masked_movie(fig, camera_paths(i).projections, mask, colors, [], ...
          {'x', 'LineWidth', 2}, {'o', 'LineWidth', 1}, 10);

      fprintf('Any key to continue...\n');
      pause;
      if ishandle(fig)
        close(fig);
      end
    end

    % Build linear systems of algebraic error.
    fprintf('Building linear systems...\n');
    proj_eqns = tracks_to_equations(cameras, tracks);

    fprintf('Solving for structure...\n');
    solution = find_structure_centroid(proj_eqns, true, settings);
    solution = shiftdim(reshape(solution, [3, num_frames, num_tracks]), 1);
    solutions(i, j).points = solution;

    if interactive
      % Plot projections with missing data.
      fig = figure();
      plot_masked_movie(fig, solution, mask, colors, [], ...
          {'x', 'LineWidth', 2}, {'o', 'LineWidth', 1}, 10);

      fprintf('Any key to continue...\n');
      pause;
      if ishandle(fig)
        close(fig);
      end
    end
  end
end
