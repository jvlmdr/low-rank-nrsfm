% S = structure_to_matrix(structure)
%
% structure -- 3 x P x F
% S -- 3F x P

function S = structure_to_matrix(structure)
  P = size(structure, 2);
  F = size(structure, 3);

  % [3, P, F] -> [3, F, P] -> [3F, P]
  S = reshape(permute(structure, [1, 3, 2]), [3 * F, P]);
end
