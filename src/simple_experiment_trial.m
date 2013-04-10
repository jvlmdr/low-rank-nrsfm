function solutions = simple_experiment_trial(solvers, scene, ranks)
  % solvers
  % - find_cameras
  % - find_cameras_with_K
  % - find_structure
  % - find_structure_with_K

  num_ranks = numel(ranks);
  num_cameras_sans_K = numel(solvers.find_cameras);
  num_cameras_with_K = numel(solvers.find_cameras_with_K);
  num_structure_sans_K = numel(solvers.find_structure);
  num_structure_with_K = numel(solvers.find_structure_with_K);

  num_cameras = num_cameras_sans_K + num_cameras_with_K;
  num_structure = num_structure_sans_K + num_structure_with_K;

  projections = scene.projections;
  % All methods get with centered projections.
  projections = bsxfun(@minus, projections, mean(projections, 2));

  % First solve for camera using methods which don't depend on K.
  for i = 1:num_cameras_sans_K
    camera_solver = solvers.find_cameras(i);
    fprintf('\nSolving for cameras using ''%s''\n', camera_solver.name);
    % Do it!
    cameras = camera_solver.solve(projections);

    for j = 1:num_structure_sans_K
      structure_solver = solvers.find_structure(j);
      fprintf('\nSolving for structure using ''%s''\n', structure_solver.name);
      fprintf('Cameras estimated using ''%s''\n', camera_solver.name);
      % Do it!
      structure = structure_solver.solve(projections, cameras);
      solution = struct('rotations', cameras, 'structure', structure);
      solutions_sans_sans(i, j) = solution;
    end

    for k = 1:num_ranks
      K = ranks(k);
      for j = 1:num_structure_with_K
        structure_solver = solvers.find_structure_with_K(j);
        fprintf('\nSolving for structure using ''%s'' with K = %d\n', ...
            structure_solver.name, K);
        fprintf('Cameras estimated using ''%s''\n', camera_solver.name);
        % Do it!
        structure = structure_solver.solve(projections, cameras, K);
        solution = struct('rotations', cameras, 'structure', structure);
        solutions_sans_with(k, i, j) = solution;
      end
    end
  end

  % Then solve for cameras using methods which do depend on K.
  for k = 1:num_ranks
    K = ranks(k);
    for i = 1:num_cameras_with_K
      camera_solver = solvers.find_cameras_with_K(i);
      fprintf('\nSolving for cameras using ''%s'' with K = %d\n', ...
          camera_solver.name, K);
      fprintf('Cameras estimated using ''%s''\n', camera_solver.name);
      % Do it!
      cameras = camera_solver.solve(projections, K);

      for j = 1:num_structure_sans_K
        structure_solver = solvers.find_structure(j);
        fprintf('\nSolving for structure using ''%s''\n', ...
            structure_solver.name);
        fprintf('Cameras estimated using ''%s''\n', camera_solver.name);
        % Do it!
        structure = structure_solver.solve(projections, cameras);
        solution = struct('rotations', cameras, 'structure', structure);
        solutions_with_sans(k, i, j) = solution;
      end

      for j = 1:num_structure_with_K
        structure_solver = solvers.find_structure_with_K(j);
        fprintf('\nSolving for structure using ''%s'' with K = %d\n', ...
            structure_solver.name, K);
        fprintf('Cameras estimated using ''%s''\n', camera_solver.name);
        % Do it!
        structure = structure_solver.solve(projections, cameras, K);
        solution = struct('rotations', cameras, 'structure', structure);
        solutions_with_with(k, i, j) = solution;
      end
    end
  end

  solutions = struct(...
      'sans_sans', solutions_sans_sans, ...
      'sans_with', solutions_sans_with, ...
      'with_sans', solutions_with_sans, ...
      'with_with', solutions_with_with);
end
