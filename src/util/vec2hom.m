function Y = vec2hom(X, dim)

if nargin < 2
  dim = 1;
end

X = shiftdim(X, dim - 1);
sz = size(X);
sz(1) = 1;
Y = cat(1, X, ones(sz));
Y = shiftdim(Y, length(sz) - (dim - 1));

end
