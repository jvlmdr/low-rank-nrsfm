% Parameters:
% X -- 3 x n

function y = nonplanarity(X)
  X = bsxfun(@minus, X, mean(X, 2));
  y = nonsubspaceness(X);
end
