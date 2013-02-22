% Parameters:
% P -- 4 x 3 projection matrix.
% X -- 3 x N matrix of N 3D points.
%
% Returns:
% 2 x N matrix of N 2D points.

function W = project(P, X)
  X = vec2hom(X);
  W = P * X;
  W = hom2vec(W);
end
