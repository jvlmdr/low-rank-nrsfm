function equations = projections_to_equations(tracks, cameras)
  num_frames = size(cameras, 3);
  num_points = length(tracks);

  % Convert from cameras and projections to equations.
  equations = struct('num_frames', num_frames);

  for i = 1:num_points
    projection_track = tracks(i);
    m = length(projection_track.frames);
    equation_track = struct('frames', projection_track.frames);

    for j = 1:m
      t = projection_track.frames(j);
      P = cameras(:, :, t);
      w = projection_track.points(j, :)';

      [A, b] = projection_to_equation(P, w);

      equation_track.equations(j) = struct('A', A, 'b', b);
      equations.tracks(i) = equation_track;
    end
  end
end
