function A = first_difference_matrix(n)
  A = toeplitz([-1; sparse(n - 2, 1)], [-1; 1; sparse(n - 2, 1)]);
end
