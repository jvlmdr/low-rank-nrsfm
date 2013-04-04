function [names, ids] = nrsfm_solver_names(solvers)
  num_full_solvers = length(solvers.nrsfm_solvers);
  num_camera_solvers = length(solvers.camera_solvers);
  num_solvers_given_cameras = length(solvers.nrsfm_solvers_given_cameras);
  num_structure_solvers = length(solvers.full_init_solvers);
  num_refiners = length(solvers.nrsfm_solvers_full_init);

  names = {};
  ids = {};

  % First solve all methods which don't require anything.
  for i = 1:num_full_solvers
    nrsfm_solver = solvers.nrsfm_solvers(i);
    names = [names, {nrsfm_solver.name}];
    ids = [ids, {nrsfm_solver.id}];
  end

  for i = 1:num_camera_solvers
    camera_solver = solvers.camera_solvers(i);

    for j = 1:num_solvers_given_cameras
      nrsfm_solver = solvers.nrsfm_solvers_given_cameras(j);
      name = sprintf('%s [%s]', nrsfm_solver.name, camera_solver.name);
      id = sprintf('%s-%s', nrsfm_solver.id, camera_solver.id);
      names = [names, {name}];
      ids = [ids, {id}];
    end

    for j = 1:num_structure_solvers
      structure_solver = solvers.full_init_solvers(j);

      for k = 1:num_refiners
        nrsfm_solver = solvers.nrsfm_solvers_full_init(k);
        name = sprintf('%s [%s] [%s]', nrsfm_solver.name, ...
            structure_solver.name, camera_solver.name);
        id = sprintf('%s-%s-%s', nrsfm_solver.id, structure_solver.id, ...
            camera_solver.id);
        names = [names, {name}];
        ids = [ids, {id}];
      end
    end
  end
end
