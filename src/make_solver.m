% Parameters:
% solve -- Function mapping F x N x 2 matrix to F x N x 3 matrix.
% name -- Nice name for printing.
% id -- Name for saving, etc.

function solver = make_solver(solve, name, id)
  if nargin == 0
    solve = [];
    name = '';
    id = '';
  end

  solver = struct('solve', solve, 'name', name, 'id', id);
end
