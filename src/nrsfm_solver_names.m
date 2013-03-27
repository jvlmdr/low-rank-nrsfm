function [names, ids] = nrsfm_solver_names(solvers)
  names = { solvers.nrsfm_solvers.name };
  ids = { solvers.nrsfm_solvers.id };

  for i = 1:length(solvers.camera_solvers)
    camera_solver = solvers.camera_solvers(i);

    camera_solver_names = cellfun(...
        @(name) { sprintf('%s [%s]', name, camera_solver.name) }, ...
        { solvers.nrsfm_solvers_given_cameras.name });
    camera_solver_ids = cellfun(...
        @(id) { sprintf('%s-%s', camera_solver.id, id) }, ...
        { solvers.nrsfm_solvers_given_cameras.id });

    names = cat(2, names, camera_solver_names);
    ids = cat(2, ids, camera_solver_ids);
  end
end
