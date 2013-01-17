% Finds S that minimizes
%   || reshape(S) ||_{*}
% subject to
%   W = R S
%
% Only applicable for parallel projection with no missing data.

function X = find_structure_matrix_equality(W, R, use_3P, settings)

  F = size(W, 1) / 2;
  N = size(W, 2);

  % Note: There might be some gain in treating this as a special case.
  % (Only having to compute the decomposition of R once.)
  % Never mind for now.

  projections = struct('Q', {}, 'u', {});
  for i = 1:N
    projections(i).Q = R;
    projections(i).q = W(:, i);
  end

  X = find_structure(projections, use_3P, settings);
end
