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
r = max(sqrt(sum(S .* S, 2)));
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

fprintf('Generating random cameras...\n');
% Generate cameras.
extrinsics = random_cameras_on_sphere(num_frames, 3 * r, shift);
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

fprintf('Building linear system of equations...\n');
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

fprintf('Solving for structure...\n');

% ADMM settings.
settings.rho = 1;
settings.mu = 10;
settings.tau_incr = 2;
settings.tau_decr = 2;
settings.max_iter = 200;
settings.epsilon_abs = 1e-3;
settings.epsilon_rel = 1e-3;

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

fprintf('Simulating missing data...\n');
% Now simulate missing data.
fraction_missing = 0.5;
num_missing = round(fraction_missing * num_frames * num_points);
subset = randperm(num_frames * num_points);
subset = subset(1:num_missing);
mask = ones(num_frames * num_points, 1);
mask(subset) = 0;
mask = reshape(mask, [num_frames, num_points]);

if interactive
  % Plot projections with missing data.
  fig = figure();
  plot_masked_movie(fig, projections, mask, 30, [])
end

% Convert from points and mask to trajectories.
tracks = struct('frames', {}, 'points', {});
for i = 1:num_points
  frames = find(mask(:, i) ~= 0);
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

% Plot 3D structure.
fig = figure();
plot_movie(fig, X, 30, [], {});
