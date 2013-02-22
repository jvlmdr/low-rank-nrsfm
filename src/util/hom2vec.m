function Y = hom2vec(X, dim)

if nargin < 2
  dim = 1;
end

X = shiftdim(X, dim - 1);
k = 1 ./ X(end, :);
Y = X(1:end-1, :) * diag(k);

sz = size(X);
sz(1) = sz(1) - 1;
Y = reshape(Y, sz);
Y = shiftdim(Y, length(sz) - (dim - 1));

end
