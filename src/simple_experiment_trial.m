function solutions = experiment_trial(solvers, scene)
  projections = scene.projections;
  projections = bsxfun(@minus, projections, mean(projections, 2));

  num_full_solvers = length(solvers.nrsfm_solvers);
  num_camera_solvers = length(solvers.camera_solvers);
  num_solvers_given_cameras = length(solvers.nrsfm_solvers_given_cameras);
  num_structure_solvers = length(solvers.full_init_solvers);
  num_refiners = length(solvers.nrsfm_solvers_full_init);

  %num_solvers = num_full_solvers + ...
  %    num_camera_solvers * (num_solvers_given_cameras + ...
  %      num_structure_solvers * num_refiners);

  solutions = [];

  % First solve all methods which don't require anything.
  for i = 1:num_full_solvers
    nrsfm_solver = solvers.nrsfm_solvers(i);
    fprintf('\nSolving ''%s''\n', nrsfm_solver.name);
    [structure, cameras] = nrsfm_solver.solve(projections);
    solution = struct('structure', structure, 'rotations', cameras);
    solutions = [solutions, solution];
  end

  for i = 1:num_camera_solvers
    camera_solver = solvers.camera_solvers(i);
    fprintf('\nSolving ''%s''\n', camera_solver.name);
    cameras_init = camera_solver.solve(projections);

    for j = 1:num_solvers_given_cameras
      nrsfm_solver = solvers.nrsfm_solvers_given_cameras(j);
      fprintf('\nSolving ''%s''\n', nrsfm_solver.name);
      fprintf('Cameras initialized by ''%s''\n', camera_solver.name);
      [structure, cameras] = nrsfm_solver.solve(projections, cameras_init);
      solution = struct('structure', structure, 'rotations', cameras);
      solutions = [solutions, solution];
    end

    for j = 1:num_structure_solvers
      structure_solver = solvers.full_init_solvers(j);
      fprintf('\nSolving ''%s''\n', structure_solver.name);
      fprintf('Cameras initialized by ''%s''\n', camera_solver.name);
      [structure_init, cameras_reinit, basis_init, coeff_init] = ...
          structure_solver.solve(projections, cameras_init);

      for k = 1:num_refiners
        nrsfm_solver = solvers.nrsfm_solvers_full_init(k);
        fprintf('\nSolving ''%s''\n', nrsfm_solver.name);
        fprintf('Cameras initialized by ''%s''\n', camera_solver.name);
        fprintf('Structure initialized by ''%s''\n', structure_solver.name);
        [structure, cameras] = nrsfm_solver.solve(projections, ...
            structure_init, cameras_reinit, basis_init, coeff_init);
        solution = struct('structure', structure, 'rotations', cameras);
        solutions = [solutions, solution];
      end
    end
  end
end
