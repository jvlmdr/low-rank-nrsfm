function solutions = noise_experiment_trial(solvers, scene)
  projections = scene.projections;
  projections = bsxfun(@minus, projections, mean(projections, 2));

  % Solve for cameras using all methods.
  cameras = arrayfun(@(solver) find_cameras(solver, projections), ...
      solvers.camera_solvers);

  % Solve different methods of full initialization given cameras.
  full_inits = arrayfun(...
      @(camera) { arrayfun(...
        @(solver) find_full_init(solver, projections, camera), ...
        solvers.full_init_solvers) }, ...
      cameras);
  full_inits = cell2mat(full_inits);

  % Append (solvers with initialization x methods of initialization) to list.
  nrsfm_solvers_given_cameras = arrayfun(...
      @(camera) { arrayfun(...
        @(solver) nrsfm_solver_given_camera(solver, camera), ...
        solvers.nrsfm_solvers_given_cameras) }, ...
      cameras);
  nrsfm_solvers_given_cameras = cell2mat(nrsfm_solvers_given_cameras);

  % Append (solvers with initialization x methods of initialization) to list.
  nrsfm_solvers_full_init = arrayfun(...
      @(full_init) { arrayfun(...
        @(solver) nrsfm_solver_full_init(solver, full_init), ...
        solvers.nrsfm_solvers_full_init) }, ...
      full_inits);
  nrsfm_solvers_full_init = cell2mat(nrsfm_solvers_full_init);

  % Concatenate.
  nrsfm_solvers = [solvers.nrsfm_solvers, nrsfm_solvers_given_cameras, ...
      nrsfm_solvers_full_init];

  % Solve!
  solutions = arrayfun(@(solver) nrsfm(solver, projections), nrsfm_solvers);
end

function solution = find_cameras(solver, projections)
  fprintf('Solving ''%s''...\n', solver.name);
  rotations = solver.solve(projections);
  solution = struct('rotations', rotations, 'solver', solver);
end

function solution = find_full_init(solver, projections, camera)
  fprintf('Solving ''%s''...\n', solver.name);
  rotations = camera.rotations;
  [structure, rotations, basis, coeff] = solver.solve(projections, rotations);
  solution = struct('rotations', rotations, 'structure', structure, ...
      'basis', basis, 'coeff', coeff, 'solver', solver);
end

function solver = nrsfm_solver_given_camera(solver, camera)
  rotations = camera.rotations;
  solve = @(projections) solver.solve(projections, rotations);
  name = [solver.name, ' [', camera.solver.name, ']'];
  id = [solver.id, '-', camera.solver.id];
  solver = make_solver(solve, name, id);
end

function solver = nrsfm_solver_full_init(solver, full_init)
  structure = full_init.structure;
  rotations = full_init.rotations;
  basis = full_init.basis;
  coeff = full_init.coeff;

  solve = @(projections) solver.solve(projections, structure, rotations, ...
      basis, coeff);
  name = [solver.name, ' [', full_init.solver.name, ']'];
  id = [solver.id, '-', full_init.solver.id];

  solver = make_solver(solve, name, id);
end

function solution = nrsfm(solver, projections)
  fprintf('Solving ''%s''...\n', solver.name);
  [structure, rotations] = solver.solve(projections);
  solution = struct('structure', structure, 'rotations', rotations);
end
