% [basis, coeff] = normalize_basis(basis, coeff)
%
% basis -- 3 x K x P
% coeff -- K x F

function [basis, coeff] = normalize_basis(basis, coeff)
  % Normalize the K basis vectors.
  b = sqrt(sum(sum(basis .^ 2, 1), 3));
  basis = bsxfun(@rdivide, basis, b);
  coeff = bsxfun(@times, coeff, b(:));
end
