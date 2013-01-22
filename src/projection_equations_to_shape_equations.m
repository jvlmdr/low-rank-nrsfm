function systems = projection_equations_to_shape_equations(projections)
  N = length(projections.tracks);
  F = projections.num_frames;

  % Matlab incredibly slow at accessing nested structs in nested for loop.
  % Extract from struct using cells and blkdiag.
  systems = struct('A', {}, 'b', {});

  for t = 1:F
    A = sparse(0, 3 * N);
    b = [];
    systems(t) = struct('A', A, 'b', b);
  end

  for i = 1:N
    track = projections.tracks(i);
    m = length(track.frames);

    for j = 1:m
      t = track.frames(j);
      Q = track.equations(j).A;
      u = track.equations(j).b;

      A = sparse(2, 3 * N);
      A(:, 3 * (i - 1) + (1:3)) = Q;
      b = u;

      systems(t).A = [systems(t).A; A];
      systems(t).b = [systems(t).b; b];
    end
  end
end
