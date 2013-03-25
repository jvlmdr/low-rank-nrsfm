% basis = basis_from_matrix(B)
%
% B -- 3K x P
% basis -- 3 x K x P

function basis = basis_from_matrix(B)
  K = size(B, 1) / 3;
  P = size(B, 2);

  % [3K, P] -> [3, K, P]
  basis = reshape(B, [3, K, P]);
end
