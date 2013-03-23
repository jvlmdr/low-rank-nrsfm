% Parameters:
% projections -- 2 x P x F
% structure -- 3 x P x F
% rotations -- 2 x 3 x F
%
% Returns:
% rotations -- 2 x 3 x F

function rotations = find_non_worse_cameras(projections, structure, rotations)
  F = size(projections, 3);

  for t = 1:F
    W_t = projections(:, :, t);
    S_t = structure(:, :, t);

    Q_old = rotations(:, :, t);
    err_old = norm(W_t - Q_old * S_t, 'fro');

    % minimize || W_t - Q S_t || s.t. Q' Q = I
    Q = procrustes(S_t', W_t')';
    err_new = norm(W_t - Q * S_t, 'fro');

    % Protect against numerical error.
    if err_new > err_old
      Q = Q_old;
    end

    rotations(:, :, t) = Q;
  end
end
