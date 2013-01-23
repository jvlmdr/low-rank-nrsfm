% Convert from cameras and projections to equations.
equations = projections_to_equations(tracks, cameras);

num_solvers = length(solvers);
solutions = struct('points', {});

for i = 1:num_solvers
  solver = solvers(i);

  fprintf('Solving for structure...\n');
  solution = solver.f(equations);
  solution = shiftdim(reshape(solution, [3, num_frames, num_points]), 1);

  solutions(i).points = solution;
end
