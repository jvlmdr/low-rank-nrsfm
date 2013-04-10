function cameras = unstack_cameras(R)
  F = size(R, 1) / 2;
  % [2F, 3] -> [2, F, 3] -> [2, 3, F]
  cameras = permute(reshape(R, [2, F, 3]), [1, 3, 2]);
end
