function solutions = noise_experiment_trial(solvers, scene)
  camera_solvers = solvers.camera_solvers;
  projections = scene.projections;

  % Solve for cameras using all methods.
  cameras = arrayfun(@(solver) find_cameras(solver, projections), ...
      camera_solvers);

  nrsfm_solvers_needing_cameras = solvers.nrsfm_solvers_given_cameras;
  % Append (solvers with initialization x methods of initialization) to list.
  nrsfm_solvers_given_cameras = arrayfun(...
      @(camera) { arrayfun(...
        @(solver) nrsfm_solver_given_camera(solver, camera), ...
        nrsfm_solvers_needing_cameras) }, ...
      cameras);
  % Take out of cells.
  nrsfm_solvers_given_cameras = cell2mat(nrsfm_solvers_given_cameras);

  % Concatenate.
  nrsfm_solvers = [solvers.nrsfm_solvers, nrsfm_solvers_given_cameras];

  % Solve!
  solutions = arrayfun(@(solver) nrsfm(solver, projections), nrsfm_solvers);
end

function solution = find_cameras(solver, projections)
  fprintf('Solving ''%s''...\n', solver.name);
  rotations = solver.solve(projections);
  solution = struct('rotations', rotations, 'solver', solver);
end

function solver = nrsfm_solver_given_camera(solver, camera)
  rotations = camera.rotations;
  solve = @(projections) solver.solve(projections, rotations);
  name = [solver.name, ' [', camera.solver.name, ']'];
  id = [solver.id, '-', camera.solver.id];
  solver = make_solver(solve, name, id);
end

function solution = nrsfm(solver, projections)
  fprintf('Solving ''%s''...\n', solver.name);
  [structure, rotations] = solver.solve(projections);
  solution = struct('structure', structure, 'rotations', rotations);
end
