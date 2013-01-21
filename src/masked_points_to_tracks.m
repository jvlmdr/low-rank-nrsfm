% Parameters:
% points -- num_frames x num_points x num_dims (usually 2 or 3)
% mask -- num_frames x num_points

function tracks = masked_points_to_tracks(points, mask);
  num_points = size(points, 2);
  num_dims = size(points, 3);

  % Convert from points and mask to tracks.
  tracks = struct('frames', {}, 'points', {});
  for i = 1:num_points
    frames = find(mask(:, i));
    m = length(frames);
    point_subset = reshape(points(frames, i, :), [m, num_dims]);

    tracks(i).frames = frames;
    tracks(i).points = point_subset;
  end
end
