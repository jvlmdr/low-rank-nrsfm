function e = projection_error(projections, structure, cameras)
  W = projections_to_matrix(projections);
  R = block_diagonal_cameras(cameras);
  S = structure_to_matrix(structure);
  e = norm(W - R * S, 'fro');
end
