function equations = tracks_to_equations(cameras, tracks)
  num_tracks = length(tracks);

  % Build linear systems of algebraic error.
  equations = struct('Q', {}, 'q', {});

  for i = 1:num_tracks
    [Q, q] = track_to_equations(cameras, tracks(i));
    equations(i).Q = Q;
    equations(i).q = q;
  end
end
