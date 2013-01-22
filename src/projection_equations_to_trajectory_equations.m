function systems = projection_equations_to_trajectory_equations(projections)
  N = length(projections.tracks);
  F = projections.num_frames;

  % Matlab incredibly slow at accessing nested structs in nested for loop.
  % Extract from struct using cells and blkdiag.
  systems = struct('A', {}, 'b', {});

  for i = 1:N
    track = projections.tracks(i);
    m = length(track.frames);

    % Build block-diagonal matrix.
    A = { track.equations.A };
    A = cellfun(@(X) { sparse(X) }, A);
    A = blkdiag(A{:});

    % Join adjacent column triples.
    Q = reshape(A, 6 * m, m);
    % Pad to size F.
    A = sparse(6 * m, F);
    A(:, track.frames) = Q;
    A = reshape(A, [2 * m, 3 * F]);

    b = { track.equations.b };
    b = cell2mat(b');

    systems(i).A = A;
    systems(i).b = b;
  end
end
