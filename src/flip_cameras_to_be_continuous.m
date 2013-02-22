% Parameters:
% R -- 2 x 3 x F

function R = flip_cameras_to_be_continuous(R)
  F = size(R, 3);

  % Negate where required.
  for t = 1:F - 1
    k1 = cross(R(1, :, t), R(2, :, t));
    k2 = cross(R(1, :, t + 1), R(2, :, t + 1));

    if dot(k1, k2) < 0
      warning('Flipping sign of cameras');
      R(:, :, t + 1) = -R(:, :, t + 1);
    end
  end
end
