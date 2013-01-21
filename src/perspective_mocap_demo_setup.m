%input_file = '../data/akhter-2008/stretch.mat';
input_file = '../data/cmu-mocap/09_09.mat';

movie_dir = '../visualize/09-09';
image_size = [848, 480];
dpi = 100;

if ~exist('interactive', 'var')
  fprintf('Defaulting to non-interactive mode\n');
  interactive = false;
end

% Load 3D points.
fprintf('Loading mocap data...\n');
%points = load_akhter_mocap(input_file);
points = load_mocap(input_file);
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

% Simulate missing data.
fprintf('Simulating missing data...\n');
random_mask = random_missing_data(num_frames, num_points, 0.5);
occlusion_mask = simulate_occlusion(num_frames, num_points, 0.5, 64);

[correspondence_mask, correspondence_index] = separate_occlusions(...
    occlusion_mask);
num_occluded_points = length(correspondence_index);

kinds_of_missing_data = [...
    struct(...
      'index', 1:num_points, ...
      'mask', ones(num_frames, num_points), ...
      'id', 'none'), ...
    struct(...
      'index', 1:num_points, ...
      'mask', random_mask, ...
      'id', 'random'), ...
    struct(...
      'index', 1:num_points, ...
      'mask', occlusion_mask, ...
      'id', 'occlusion'), ...
    struct(...
      'index', correspondence_index, ...
      'mask', correspondence_mask, ...
      'id', 'correspondence')];

% Plot sparsity patterns.
fig = figure();
subplot(2, 2, 1);
imshow(1 - random_mask);
title('Random missing data');
subplot(2, 2, 2);
imshow(1 - occlusion_mask);
title('Simulated occlusion');
subplot(2, 2, 3);
imshow(1 - correspondence_mask);
title('Correspondence loss');

fprintf('Generating random cameras...\n');
smooth_extrinsics = smooth_random_cameras_on_sphere(num_frames, 3 * radius, ...
    5 * pi / 180, radius, 1 * pi / 180, shift);
random_extrinsics = random_cameras_on_sphere(num_frames, 3 * radius, shift);

% Generate cameras.
camera_paths = [...
    struct(...
      'extrinsics', smooth_extrinsics, ...
      'id', 'smooth'), ...
    struct(...
      'extrinsics', random_extrinsics, ...
      'id', 'random')];

for i = 1:length(camera_paths)
  extrinsics = camera_paths(i).extrinsics;

  K = diag([1, 1, -1]);
  cameras = zeros(3, 4, num_frames);
  for t = 1:num_frames
    cameras(:, :, t) = K * extrinsics(:, :, t);
  end

  camera_paths(i).cameras = cameras;
  camera_paths(i).projections = project(cameras, points);
end

% Plot cameras.
fig = figure();
title('Cameras');
subplot(1, 2, 1);
plot_cameras(fig, camera_paths(1).extrinsics, 1);
subplot(1, 2, 2);
plot_cameras(fig, camera_paths(2).extrinsics, 1);

% Generate projections for each configuration.
sequences = struct('tracks', {});

for i = 1:length(camera_paths)
  % Project points on to cameras.
  cameras = camera_paths(i).cameras;
  projections = camera_paths(i).projections;

  for j = 1:length(kinds_of_missing_data)
    index = kinds_of_missing_data(j).index;
    mask = kinds_of_missing_data(j).mask;

    sequences(i, j) = struct(...
        'tracks', masked_points_to_tracks(projections(:, index, :), mask));
  end
end

colors = hsv(num_occluded_points);
colors = colors(randperm(num_occluded_points), :);
