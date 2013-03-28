function x = nuclear_norm(X)
  x = sum(svd(X));
end
