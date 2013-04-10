function R = stack_cameras(cameras)
  F = size(cameras, 3);
  % [2, 3, F] -> [2, F, 3] -> [2F, 3]
  R = reshape(permute(cameras, [1, 3, 2]), [2 * F, 3]);
end
