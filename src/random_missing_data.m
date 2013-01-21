% Parameters:
% fraction -- Fraction of points missing, between 0 and 1.

function mask = random_missing_data(num_frames, num_points, fraction)
  n = num_frames * num_points;
  % Number of points missing.
  num_visible = round((1 - fraction) * n);
  subset = randperm(n);
  subset = subset(1:num_visible);
  mask = sparse(n, 1);
  mask(subset) = 1;
  mask = reshape(mask, [num_frames, num_points]);
end
