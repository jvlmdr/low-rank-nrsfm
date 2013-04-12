function solutions = exact_rank_experiment_trial(solvers, scene)
  % solvers
  % - find_both
  % - find_cameras
  % - find_structure

  num_both = numel(solvers.find_both);
  num_cameras = numel(solvers.find_cameras);
  num_structure = numel(solvers.find_structure);

  projections = scene.projections;
  % All methods get with centered projections.
  projections = bsxfun(@minus, projections, mean(projections, 2));

  F = size(projections, 3);

  % Extract cameras, remove scales.
  true_cameras = zeros(2, 3, F);
  for t = 1:F
    R = scene.cameras(t).P(1:2, 1:3);
    scale = norm(R, 'fro') / sqrt(2);
    true_cameras(:, :, t) = 1 / scale * R;
  end

  solutions_joint(num_both) = ...
      struct('rotations', [], 'structure', []);
  solutions_independent(num_cameras, num_structure) = ...
      struct('rotations', [], 'structure', []);

  % Solve methods which estimate everything from projections.
  for i = 1:num_both
    solver = solvers.find_both(i);
    fprintf('\nSolving for structure and cameras using ''%s''\n', solver.name);
    % Do it!
    [structure, cameras] = solver.solve(projections);

    solution = struct('rotations', cameras, 'structure', structure);
    solutions_joint(i) = solution;
  end

  for i = 1:num_cameras
    camera_solver = solvers.find_cameras(i);
    fprintf('\nSolving for cameras using ''%s''\n', camera_solver.name);
    % Do it!
    cameras = camera_solver.solve(projections, true_cameras);

    % Even camera solvers have a bad day sometimes...
    if ~isempty(cameras)
      for j = 1:num_structure
        structure_solver = solvers.find_structure(j);
        fprintf('\nSolving for structure using ''%s''\n', ...
            structure_solver.name);
        fprintf('Cameras estimated using ''%s''\n', camera_solver.name);
        % Do it!
        structure = structure_solver.solve(projections, cameras);

        solution = struct('rotations', cameras, 'structure', structure);
        solutions_independent(i, j) = solution;
      end
    end
  end

  solutions = struct(...
      'joint', solutions_joint, ...
      'independent', solutions_independent);
end
