% structure = compose_structure(basis, coeff)
%
% Parameters:
% basis -- 3 x K x P
% coeff -- K x F
%
% Returns:
% structure -- 3 x P x F

function structure = compose_structure(basis, coeff)
  K = size(basis, 2);
  P = size(basis, 3);
  F = size(coeff, 2);

  % [3, K, P] -> [3, P, K] -> [3P, K]
  B = reshape(permute(basis, [1, 3, 2]), [3 * P, K]);
  C = coeff;

  S = B * C;

  % [3P, F] -> [3, P, F]
  structure = reshape(S, [3, P, F]);
end
