function structure = find_planar_structure(projections, cameras)
  P = size(projections, 2);
  F = size(projections, 3);

  structure = zeros(3, P, F);

  for t = 1:F
    structure(:, :, t) = pinv(cameras(:, :, t)) * projections(:, :, t);
  end
end
