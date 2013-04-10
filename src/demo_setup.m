close all;
rng('default');

K_range = [2, 4, 6];
scale_stddev = 1; %sqrt(2);
omega_stddev = 5 / 180 * pi;
use_akhter_data = input('Use data of Akhter? (true, false) ');

if use_akhter_data
  % Load mocap sequence.
  seq_name = input('Sequence? (drink, pickup, stretch, yoga) ', 's');
  mocap_file = ['../data/akhter-2008/', seq_name];
  load(mocap_file, 'S', 'Rs');
  F = size(Rs, 1) / 2;
  P = size(S, 2);

%  % [2F, 3] -> [2, F, 3] -> [2, 3, F]
%  rotations = permute(reshape(Rs, [2, F, 3]), [1, 3, 2]);
%  % Scale each frame.
%  scales = exp(log(scale_stddev) * randn(F, 1));
%  cameras = bsxfun(@times, rotations, reshape(scales, [1, 1, F]));
  clear Rs;

  world_structure = structure_from_matrix(S);
  clear S;
else
  % Load mocap sequence.
  data = load('../data/mocap-data.mat');
  num_sequences = size(data.sequences, 4);
  mocap_index = input(sprintf('Mocap index? (1, ..., %d) ', num_sequences));
  F = size(data.sequences, 1);
  P = size(data.sequences, 2);
  world_structure = data.sequences(:, :, :, mocap_index);
  % [F, P, 3] -> [3, P, F]
  world_structure = permute(world_structure, [3, 2, 1]);
end

% Angular change in each frame.
%omegas = omega_stddev * randn(F, 1);
omegas = omega_stddev * ones(F, 1);
% Angle in each frame.
thetas = cumsum(omegas) + rand() * 2 * pi;
% Scale in each frame.
scales = exp(log(scale_stddev) * randn(F, 1));

% Generate camera motion.
scene = generate_scene_for_sequence(world_structure, thetas, scales);

% Extract cameras, remove scales.
scales = zeros(F, 1);
world_cameras = zeros(2, 3, F);
for t = 1:F
  R = scene.cameras(t).P(1:2, 1:3);
  scales(t) = norm(R, 'fro') / sqrt(2);
  world_cameras(:, :, t) = 1 / scales(t) * R;
end

% Subtract centroid.
centroid = mean(world_structure, 2);
unscaled_unaligned_structure = bsxfun(@minus, world_structure, centroid);
% Apply scale to structure instead of cameras.
unaligned_structure = bsxfun(@times, unscaled_unaligned_structure, ...
    reshape(scales, [1, 1, F]));

% Project.
R = block_diagonal_cameras(world_cameras);
S = structure_to_matrix(unaligned_structure);
W = R * S;
projections = projections_from_matrix(W);

K = input('Choose K? ');

% Align shapes.
align = input('Align using congealing? ');
if align
  unscaled_structure = congeal_shapes(unscaled_unaligned_structure, 1e-6, 20);
  %unscaled_structure = congeal_shapes_low_rank(unscaled_structure, K, 1e-6, 20);
  % Calculate rotations after the fact.
  rotations = zeros(3, 3, F);
  for t = 1:F
    A = unscaled_unaligned_structure(:, :, t);
    B = unscaled_structure(:, :, t);
    R = procrustes(A', B')';
    rotations(:, :, t) = R;
  end
else
  unscaled_structure = unscaled_unaligned_structure;
  rotations = repmat(eye(3, 3), [1, 1, F]);
end

% Apply inverse rotations to cameras.
true_cameras = zeros(2, 3, F);
for t = 1:F
  R = rotations(:, :, t);
  true_cameras(:, :, t) = world_cameras(:, :, t) * R';
end
% Apply scale to structure.
true_structure = bsxfun(@times, unscaled_structure, reshape(scales, [1, 1, F]));

