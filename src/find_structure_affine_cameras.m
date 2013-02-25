% Finds S that minimizes
%   || reshape(S) ||_{*}
% subject to
%   w_ti = R_t s_ti.
%
% Parameters:
% W -- 2F x P
% R -- 2 x 3 x F

function X = find_structure_affine_cameras(W, R, use_3P, settings)
  N = size(W, 2);
  F = size(R, 3);

  W = reshape(W, [2, F, N]);

  % Build systems of equations.
  projections = struct('num_frames', F);
  for i = 1:N
    track = struct('frames', [], 'equations', struct([]));
    for t = 1:F
      track.frames(t) = t;
      track.equations(t).A = R(:, :, t);
      track.equations(t).b = W(:, t, i);
    end
    projections.tracks(i) = track;
  end

  X = find_structure(projections, use_3P, settings);
end
