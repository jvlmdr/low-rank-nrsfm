% projections = projections_from_matrix(W)
%
% W -- 2F x P
% projections -- 2 x P x F

function projections = projections_from_matrix(W)
  F = size(W, 1) / 2;
  P = size(W, 2);

  % [2F, P] -> [2, F, P] -> [2, P, F]
  projections = permute(reshape(W, [2, F, P]), [1, 3, 2]);
end
