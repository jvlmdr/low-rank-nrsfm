output_dir = ['../output/florence/'];

load('../data/florence/cameras');

num_frames = length(C);
num_points = size(TrackM, 2);
camera_size = 0.2;

% Load extrinsic projection matrices into tensor.
cameras = zeros(3, 4, num_frames);
for t = 1:num_frames
  cameras(:, :, t) = C{t}.P;
end

% Convert to trajectories data structure.
projections = struct('num_frames', num_frames);

for i = 1:num_points
  projections.tracks(i) = struct('frames', [], 'points', []);
end
for t = 1:num_frames
  for j = 1:size(C{t}.m, 1);
    i = C{t}.m(j, 3);
    track = projections.tracks(i);
    track.frames = [track.frames; t];
    track.points = [track.points; C{t}.m(j, 1:2)];
    projections.tracks(i) = track;
  end
end

% Plot cameras.
fig = figure();
plot_cameras(fig, cameras, 0.1);
fprintf('Press enter to continue...\n');
pause;
close(fig);

%% Render movie of cameras.
%frame_dir = [output_dir, 'cameras/'];
%unix(['rm -rf ', frame_dir]);
%unix(['mkdir -p ', frame_dir]);
%frame_format = [frame_dir, '%07d.png'];
%video = [output_dir, 'cameras.mp4'];
%fig = figure();
%save_camera_movie(fig, cameras, ...
%    @(fig) view(gca(fig), [0, 90]), ...
%    @(fig, t) print_image(fig, [848, 480], 120, ...
%      sprintf(frame_format, t), {'-dpng'}), ...
%    0.1);
%close(fig);
%unix(['ffmpeg -sameq -y -i ', frame_format, ' ', video]);

% Convert from cameras and projections to equations.
equations = projections_to_equations(projections.tracks, cameras);

% Stock standard ADMM settings.
settings.rho = 1;
settings.mu = 10;
settings.tau_incr = 2;
settings.tau_decr = 2;
settings.max_iter = 100;
settings.epsilon_abs = 1e-3;
settings.epsilon_rel = 1e-3;

solvers = [...
  struct(...
    'f', @(equations) find_structure(equations, false, settings), ...
    'name', 'Nuclear norm (3FxP)', ...
    'id', 'nuclear_3F'), ...
  struct(...
    'f', @(equations) find_structure(equations, true, settings), ...
    'name', 'Nuclear norm (3PxF)', ...
    'id', 'nuclear_3P'), ...
  struct(...
    'f', @(equations) find_structure_centroid(equations, false, ...
      settings), ...
    'name', 'Nuclear norm (3FxP) with global centroid', ...
    'id', 'nuclear_3F_centroid'), ...
  struct(...
    'f', @(equations) find_structure_centroid(equations, true, ...
      settings), ...
    'name', 'Nuclear norm (3PxF) with global centroid', ...
    'id', 'nuclear_3P_centroid'), ...
  struct(...
    'f', @(equations) find_structure_smooth_centroid(equations, false, ...
      first_difference_matrix(num_frames), 1, settings), ...
    'name', 'Nuclear norm (3FxP) with filter prior and global centroid', ...
    'id', 'nuclear_3F_smooth_centroid'), ...
  struct(...
    'f', @(equations) find_structure_smooth_centroid(equations, true, ...
      first_difference_matrix(num_frames), 1, settings), ...
    'name', 'Nuclear norm (3PxF) with filter prior and global centroid', ...
    'id', 'nuclear_3P_smooth_centroid'), ...
  struct(...
    'f', @(equations) find_structure_differences(equations, false, ...
      settings), ...
    'name', 'Nuclear norm (3FxP) of differences', ...
    'id', 'nuclear_3F_differences'), ...
  struct(...
    'f', @(equations) find_structure_differences(equations, true, ...
      settings), ...
    'name', 'Nuclear norm (3PxF) of differences', ...
    'id', 'nuclear_3P_differences'), ...
];

num_solvers = length(solvers);

% Run each solver.
clear solutions;
for i = 1:num_solvers
  fprintf('Method %u: "%s"\n', i, solvers(i).name);
  X = solvers(i).f(equations);
  solutions(i).X = X;

  % [3F, P] to [F, P, 3]
  points = shiftdim(reshape(X, [3, num_frames, num_points]), 1);
  % Generate 2D projections.
  reprojected_points = project(cameras, points);

%  % Restrict 2D projections to observed features.
%  clear reprojections;
%  for i = 1:num_points
%    frames = find(TrackM(:, i));
%    proj = reshape(reprojected_points(frames, i, :), length(frames), 2);
%    reprojections(i) = struct('frames', frames, 'points', proj);
%  end
%
%  % Check projection error.
%  e = [];
%  for i = 1:num_points
%    assert(all(projections(i).frames == reprojections(i).frames));
%    d = projections(i).points - reprojections(i).points;
%    e = [e; sqrt(sum(d .* d, 2))];
%  end
%  fprintf('Average reprojection error => %g\n', mean(e));

  % Plot movie.
  %fig = figure();
  %plot_reconstruction_movie(fig, points, TrackM, cameras, 30, [], camera_size);
  %fprintf('Press enter to continue...\n');
  %pause;
  %close(fig);
end

% Save video for each solver.
for i = 1:num_solvers
  frame_dir = [output_dir, solvers(i).id, '/'];
  unix(['rm -rf ', frame_dir]);
  unix(['mkdir -p ', frame_dir]);
  frame_format = [frame_dir, '%07d.png'];

  points = shiftdim(reshape(solutions(i).X, [3, num_frames, num_points]), 1);
  fig = figure();
  colors = hsv(num_points);
  save_reconstruction_movie(fig, points, TrackM, colors, cameras, camera_size, ...
      [], true, {'x', 'LineWidth', 2}, {'o', 'LineWidth', 1}, [848, 480], ...
      100, output_dir, solvers(i).id);
  close(fig);

  video = [output_dir, solvers(i).id, '.mp4'];
  unix(['ffmpeg -sameq -y -i ', frame_format, ' ', video]);
  unix(['rm -rf ', frame_dir]);
end
