% Parameters:
% Rs -- 2 x 3 x F
%
% Returns:
% R -- 2F x 3F, sparse

function R = block_diagonal_cameras(Rs)
  F = size(Rs, 3);

  R = mat2cell(Rs, 2, 3, ones(F, 1));
  R = cellfun(@(X) { sparse(X) }, R);
  R = blkdiag(R{:});
end
