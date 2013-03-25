% B = basis_to_matrix(basis)
%
% basis -- 3 x K x P
% B -- 3K x P

function B = basis_to_matrix(basis)
  K = size(B, 2);
  P = size(B, 3);

  % [3, K, P] -> [3K, P]
  basis = reshape(B, [3 * K, P]);
end
