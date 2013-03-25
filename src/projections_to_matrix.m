% W = projections_to_matrix(projections)
%
% projections -- 2 x P x F
% W -- 2F x P

function W = projections_to_matrix(projections)
  P = size(projections, 2);
  F = size(projections, 3);

  % [2, P, F] -> [2, F, P] -> [2F, P]
  W = reshape(permute(projections, [1, 3, 2]), [2 * F, P]);
end
