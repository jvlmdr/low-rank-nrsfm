% Returns a linear system of equations given a camera and a projection
%
% Parameters:
% P -- 3 x 4
% w -- 2 x 1
%
% Returns:
% Q -- 2 x 3
% q -- 2 x 1
%   Q x = q defines the projection constraints

function [Q, q] = projection_to_equation(P, w)
  % P [x; 1] ~ [w 1]
  % P(:, 1:3) x + P(:, 4) ~ [w; 1]
  % P(1:2, 1:3) x + P(1:2, 4) / (P(3, 1:3) x + P(3, 4)) = w
  % P(1:2, 1:3) x + P(1:2, 4) = (P(3, 1:3) x + P(3, 4)) w
  % (P(1:2, 1:3) - w P(3, 1:3)) x = P(3, 4) w - P(1:2, 4)
  % Q x = q

  Q = P(1:2, 1:3) - w * P(3, 1:3);
  q = P(3, 4) * w - P(1:2, 4);
end
