% Returns a linear system of equations given cameras and projections
%
% Parameters:
% cameras -- 3 x 4 x F
% track
%   frames -- m x 1 vector of frame indices
%   points -- m x 2 vector of observations
%
% Returns:
% Q -- 2m x 3F sparse matrix
% q -- 2m x 1 vector
%   Q x = q defines the projection constraints

function [Q, q] = track_to_equations(cameras, track)
  F = size(cameras, 3);
  m = length(track.frames);

  Q = sparse(2 * m, 3 * F);
  q = zeros(2 * m, 1);

  for i = 1:m
    w = track.points(i, :)';
    t = track.frames(i);
    P = cameras(:, :, t);

    % P [x; 1] ~ [w 1]
    % P(:, 1:3) x + P(:, 4) ~ [w; 1]
    % P(1:2, 1:3) x + P(1:2, 4) / (P(3, 1:3) x + P(3, 4)) = w
    % P(1:2, 1:3) x + P(1:2, 4) = (P(3, 1:3) x + P(3, 4)) w
    % (P(1:2, 1:3) - w P(3, 1:3)) x = P(3, 4) w - P(1:2, 4)
    % Q x = q

    Qt = P(1:2, 1:3) - w * P(3, 1:3);
    qt = P(3, 4) * w - P(1:2, 4);

    Q(2 * (i - 1) + (1:2), 3 * (t - 1) + (1:3)) = Qt;
    q(2 * (i - 1) + (1:2)) = qt;
  end

end
