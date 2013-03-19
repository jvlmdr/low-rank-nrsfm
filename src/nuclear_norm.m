function x = nuclear_norm(X)
  return sum(svd(X));
end
