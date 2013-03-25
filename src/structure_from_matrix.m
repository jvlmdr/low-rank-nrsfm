% structure = structure_from_matrix(S)
%
% S -- 3F x P
% structure -- 3 x P x F

function structure = structure_from_matrix(S)
  F = size(S, 1) / 3;
  P = size(S, 2);

  % [3F, P] -> [3, F, P] -> [3, P, F]
  structure = permute(reshape(S, [3, F, P]), [1, 3, 2]);
end
