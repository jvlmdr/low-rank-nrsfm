image_size = [848, 480];
dpi = 100;

% Load 3D points.
fprintf('Loading mocap data...\n');

%input_file = '../data/akhter-2008/yoga.mat';
%movie_dir = '../visualize/translation-experiment/yoga';
%points = load_akhter_mocap(input_file);
input_file = '../data/cmu-mocap/09_09.mat';
movie_dir = '../visualize/translation-experiment/09-09';
points = load_mocap(input_file);

if ~exist('interactive', 'var')
  fprintf('Defaulting to non-interactive mode\n');
  interactive = false;
end

num_frames = size(points, 1);
num_points = size(points, 2);

% Subtract centroid of entire sequence.
points = subtract_mean(points);
% Compute maximum radius.
radius = max(max(sqrt(sum(points .* points, 3))));

% Add random shift of entire scene (for effects of centroid).
shift = 100 * randn(3, 1);
points = reshape(points, [num_frames * num_points, 3]);
points = points + ones(num_frames * num_points, 1) * shift';
points = reshape(points, [num_frames, num_points, 3]);

% Plot input sequence.
if interactive
  fig = figure();
  plot_movie(fig, points, 30, [], {'kx', 'LineWidth', 2});

  fprintf('Any key to continue...\n');
  pause;
  if ishandle(fig)
    close(fig);
  end
end

fprintf('Generating random cameras...\n');
% Generate camera motion.
extrinsics = smooth_random_cameras_on_sphere(num_frames, 1.5 * radius, ...
    1 * pi / 180, 0.5 * radius, 1 * pi / 180, shift);
% Incorporate intrinsics.
K = diag([1, 1, -1]);
cameras = zeros(3, 4, num_frames);
for t = 1:num_frames
  cameras(:, :, t) = K * extrinsics(:, :, t);
end

% Plot cameras.
fig = figure();
title('Cameras');
plot_cameras(fig, extrinsics, 0.2);

% Generate projections for each configuration.
projections = project(cameras, points);
tracks = points_to_tracks(projections);

% Generate a unique color for each 3D point.
colors = hsv(num_points);
colors = colors(randperm(num_points), :);

% ADMM settings.
settings = struct(...
    'rho', 1, ...
    'mu', 10, ...
    'tau_incr', 2, ...
    'tau_decr', 2, ...
    'max_iter', 100, ...
    'epsilon_abs', 1e-3, ...
    'epsilon_rel', 1e-3);

gamma = 0;
centroid = reshape(mean(points, 2), [num_frames, 3]);
for t = 1:num_frames
  k = extrinsics(3, 1:3, t)';
  gamma = gamma + k' * centroid(t, :)';
end
gamma_approx = gamma / num_frames

solvers = [...
  struct(...
    'f', @(equations) find_structure_centroid(equations, true, settings), ...
    'name', 'Global centroid', ...
    'id', 'global'), ...
  struct(...
    'f', @(equations) find_structure_free_translation(equations, true, ...
      settings), ...
    'name', 'Centroid per frame', ...
    'id', 'free'), ...
  struct(...
    'f', @(equations) find_structure_constrained_translation(equations, ...
      extrinsics, 100, true, settings), ...
    'name', 'Centroid per frame (\gamma = 100)', ...
    'id', 'constrained-100'), ...
  struct(...
    'f', @(equations) find_structure_constrained_translation(equations, ...
      extrinsics, 200, true, settings), ...
    'name', 'Centroid per frame (\gamma = 200)', ...
    'id', 'constrained-200'), ...
  struct(...
    'f', @(equations) find_structure_constrained_translation(equations, ...
      extrinsics, 500, true, settings), ...
    'name', 'Centroid per frame (\gamma = 500)', ...
    'id', 'constrained-500'), ...
];
