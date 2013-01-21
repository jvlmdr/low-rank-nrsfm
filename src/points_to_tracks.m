% Parameters:
% points -- num_frames x num_points x num_dims (usually 2 or 3)

function tracks = points_to_tracks(points, mask);
  num_frames = size(points, 1);
  num_points = size(points, 2);
  num_dims = size(points, 3);

  % Convert from points and mask to tracks.
  tracks = struct('frames', {}, 'points', {});
  for i = 1:num_points
    point_subset = reshape(points(:, i, :), [num_frames, num_dims]);
    tracks(i).frames = 1:num_frames;
    tracks(i).points = point_subset;
  end
end
