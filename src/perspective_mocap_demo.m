if ~exist('interactive', 'var')
  fprintf('Defaulting to non-interactive mode\n');
  interactive = false;
end

fprintf('Loading mocap data...\n');
% Load 3D points.
load('../data/akhter-2008/yoga');

num_frames = size(S, 1) / 3;
num_points = size(S, 2);

% [3F, N] -> [F, N, 3]
S = shiftdim(reshape(S, [3, num_frames, num_points]), 1);
% Swap back y and z.
S = S(:, :, [1, 3, 2]);

shift = 100 * randn(3, 1);

% Subtract mean.
S = reshape(S, [num_frames * num_points, 3]);
mu = mean(S)';
% Subtract mean.
S = S - ones(num_frames * num_points, 1) * mu';
% Compute maximum radius.
radius = max(sqrt(sum(S .* S, 2)));
% Add shift.
S = S + ones(num_frames * num_points, 1) * shift';
% Restore dimension.
S = reshape(S, [num_frames, num_points, 3]);

if interactive
  % Plot 3D structure.
  fig = figure();
  plot_movie(fig, S, 30, [], {});
  fprintf('Press enter to continue...\n');
  pause;
  if ishandle(fig)
    close(fig);
  end
end

fprintf('Simulating random missing data...\n');
% Now simulate missing data.
fraction_missing = 0.5;
num_missing = round(fraction_missing * num_frames * num_points);
subset = randperm(num_frames * num_points);
subset = subset(1:num_missing);
random_mask = ones(num_frames * num_points, 1);
random_mask(subset) = 0;
random_mask = reshape(random_mask, [num_frames, num_points]);

if interactive
  % Plot sparsity pattern.
  fig = figure();
  spy(random_mask);

  fprintf('Press enter to continue...\n');
  pause;
  if ishandle(fig)
    close(fig);
  end
end

fprintf('Simulating realistic missing data...\n');
% Probability of transitioning visible = 0 to 1 and 1 to 0, respectively.
p_change = [0.05; 0.05];
% Initial state.
p_visible = 0.5;
occlusion_mask = zeros(num_frames, num_points);

for i = 1:num_points
  r_init = rand();
  r = rand(num_frames - 1, 1);

  % Frames t:num_frames still need to be labelled.
  t = 1;
  v = r_init < p_visible;
  visible = zeros(0, 1);

  while t <= num_frames
    % Look for the next transition event from state v.
    j = find(r(t:end) < p_change(v + 1), 1);

    if isempty(j)
      % No changes until the end.
      j = num_frames - t + 1;
    end

    % Next j frames are the same.
    visible(t:t + j - 1) = v;
    t = t + j;
    v = 1 - v;
  end

  occlusion_mask(:, i) = visible;
end
fraction_missing = 1 - nnz(occlusion_mask) / numel(occlusion_mask);
fprintf('Missing data: %.2g%%\n', fraction_missing * 100);

if interactive
  % Plot sparsity pattern.
  fig = figure();
  spy(occlusion_mask);

  fprintf('Press enter to continue...\n');
  pause;
  if ishandle(fig)
    close(fig);
  end
end

% Break occlusions into contiguous observations.
correspondence_mask = sparse(num_frames, 0);
index = [];
n = 0;
for i = 1:num_points
  % Still need to find visible segments in frames t..num_frames.
  t = 1;
  found_all = false;

  while ~found_all
    % Find the next visible frame.
    j = find(occlusion_mask(t:num_frames, i), 1);

    if isempty(j)
      found_all = true;
    else
      t = t + j - 1;

      % There is a visible segment.
      % Find the next occluded frame.
      j = find(~occlusion_mask(t:num_frames, i), 1);

      if isempty(j)
        % Go to end of sequence.
        j = num_frames - t + 1;
      end

      % Extract.
      visible = sparse(num_frames, 1);
      visible(t:t + j - 1) = 1;
      correspondence_mask = [correspondence_mask, visible];
      index = [index, i];
      n = n + 1;
      t = t + j;
    end
  end
end

if interactive
  % Plot sparsity pattern.
  fig = figure();
  spy(correspondence_mask);

  fprintf('Press enter to continue...\n');
  pause;
  if ishandle(fig)
    close(fig);
  end
end

for smooth = [false, true]
  fprintf('Generating random cameras...\n');
  % Generate cameras.
  if smooth
    extrinsics = smooth_random_cameras_on_sphere(num_frames, 3 * radius, ...
        5 * pi / 180, radius, 1 * pi / 180, shift);
  else
    extrinsics = random_cameras_on_sphere(num_frames, 3 * radius, shift);
  end
  K = diag([1, 1, -1]);
  cameras = zeros(3, 4, num_frames);
  for t = 1:num_frames
    cameras(:, :, t) = K * extrinsics(:, :, t);
  end

  if interactive
    % Plot cameras.
    fig = figure();
    plot_cameras(fig, extrinsics, 1);
    title('Cameras');
    fprintf('Press enter to continue...\n');
    pause;
    if ishandle(fig)
      close(fig);
    end
  end

  fprintf('Projecting structure into cameras...\n');
  % Project 3D structure into cameras.
  projections = project(cameras, S);

  if interactive
    % Plot projections.
    fig = figure();
    plot_movie(fig, projections, 10, [], {});
    fprintf('Press enter to continue...\n');
    pause;
    if ishandle(fig)
      close(fig);
    end
  end

  fprintf('Building linear systems...\n');
  % Build system of equations from projection constraints.
  observations = struct('Q', {}, 'q', {});
  for i = 1:num_points
    Q = {};
    q = [];

    for t = 1:num_frames
      P = cameras(:, :, t);
      w = reshape(projections(t, i, :), [2, 1]);

      Q{t} = P(1:2, 1:3) - w * P(3, 1:3);
      q(:, t) = P(3, 4) * w - P(1:2, 4);
    end

    Q = cellfun(@(X) { sparse(X) }, Q);
    observations(i).Q = blkdiag(Q{:});
    observations(i).q = q(:);
  end

  % ADMM settings.
  settings.rho = 1;
  settings.mu = 10;
  settings.tau_incr = 2;
  settings.tau_decr = 2;
  settings.max_iter = 200;
  settings.epsilon_abs = 1e-3;
  settings.epsilon_rel = 1e-3;

  fprintf('Solving for structure...\n');
  X = find_structure_centroid(observations, true, settings);
  X = shiftdim(reshape(X, [3, num_frames, num_points]), 1);

  e = X - S;
  e = reshape(e, [num_frames * num_points, 3]);
  e = sqrt(sum(e .* e, 2));
  e = mean(e)

  if interactive
    % Plot 3D structure.
    fig = figure();
    plot_movie(fig, X, 30, [], {});
    fprintf('Press enter to continue...\n');
    pause;
    if ishandle(fig)
      close(fig);
    end
  end

  if interactive
    % Plot projections with missing data.
    fig = figure();
    plot_masked_movie(fig, projections, random_mask, 30, [])
    fprintf('Press enter to continue...\n');
    pause;
    if ishandle(fig)
      close(fig);
    end
  end

  % Convert from points and mask to trajectories.
  tracks = struct('frames', {}, 'points', {});
  for i = 1:num_points
    frames = find(random_mask(:, i) ~= 0);
    m = length(frames);
    points = reshape(projections(frames, i, :), [m, 2]);

    tracks(i).frames = frames;
    tracks(i).points = points;
  end

  fprintf('Building linear systems...\n');
  % Build linear systems of algebraic error.
  projection_constraints = struct('Q', {}, 'q', {});
  for i = 1:num_points
    [Q, q] = track_to_equations(cameras, tracks(i));
    projection_constraints(i).Q = Q;
    projection_constraints(i).q = q;
  end

  fprintf('Solving for structure...\n');
  X = find_structure_centroid(projection_constraints, true, settings);
  X = shiftdim(reshape(X, [3, num_frames, num_points]), 1);

  e = X - S;
  e = reshape(e, [num_frames * num_points, 3]);
  e = sqrt(sum(e .* e, 2));
  e = mean(e)

  if interactive
    % Plot 3D structure.
    fig = figure();
    plot_masked_movie(fig, X, random_mask, 30, []);
    fprintf('Press enter to continue...\n');
    pause;
    if ishandle(fig)
      close(fig);
    end
  end

  if interactive
    % Plot projections with missing data.
    fig = figure();
    plot_masked_movie(fig, projections, occlusion_mask, 30, [])
    fprintf('Press enter to continue...\n');
    pause;
    if ishandle(fig)
      close(fig);
    end
  end

  % Convert from points and mask to trajectories.
  tracks = struct('frames', {}, 'points', {});
  for i = 1:num_points
    frames = find(occlusion_mask(:, i) ~= 0);
    m = length(frames);
    points = reshape(projections(frames, i, :), [m, 2]);

    tracks(i).frames = frames;
    tracks(i).points = points;
  end

  fprintf('Building linear systems...\n');
  % Build linear systems of algebraic error.
  projection_constraints = struct('Q', {}, 'q', {});
  for i = 1:num_points
    [Q, q] = track_to_equations(cameras, tracks(i));
    projection_constraints(i).Q = Q;
    projection_constraints(i).q = q;
  end

  X = find_structure_centroid(projection_constraints, true, settings);
  X = shiftdim(reshape(X, [3, num_frames, num_points]), 1);

  e = X - S;
  e = reshape(e, [num_frames * num_points, 3]);
  e = sqrt(sum(e .* e, 2));
  e = mean(e)

  if interactive
    % Plot 3D structure.
    fig = figure();
    plot_masked_movie(fig, X, occlusion_mask, 30, []);
    fprintf('Press enter to continue...\n');
    pause;
    if ishandle(fig)
      close(fig);
    end
  end

  if interactive
    % Plot projections with missing data.
    fig = figure();
    plot_masked_movie(fig, projections(:, index, :), correspondence_mask, ...
        30, []);
    fprintf('Press enter to continue...\n');
    pause;
    if ishandle(fig)
      close(fig);
    end
  end

  % Convert from points and mask to trajectories.
  tracks = struct('frames', {}, 'points', {});
  num_tracks = size(correspondence_mask, 2);
  for i = 1:num_tracks
    frames = find(correspondence_mask(:, i) ~= 0);
    m = length(frames);
    points = reshape(projections(frames, index(i), :), [m, 2]);

    tracks(i).frames = frames;
    tracks(i).points = points;
  end

  fprintf('Building linear systems...\n');
  % Build linear systems of algebraic error.
  projection_constraints = struct('Q', {}, 'q', {});
  for i = 1:num_tracks
    [Q, q] = track_to_equations(cameras, tracks(i));
    projection_constraints(i).Q = Q;
    projection_constraints(i).q = q;
  end

  fprintf('Solving for structure...\n');
  X = find_structure_centroid(projection_constraints, true, settings);
  X = shiftdim(reshape(X, [3, num_frames, num_tracks]), 1);

  if interactive
    % Plot 3D structure.
    fig = figure();
    plot_masked_movie(fig, X, correspondence_mask, 30, []);
    fprintf('Press enter to continue...\n');
    pause;
    if ishandle(fig)
      close(fig);
    end
  end
end
