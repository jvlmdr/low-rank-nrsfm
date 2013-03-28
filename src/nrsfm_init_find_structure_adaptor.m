% Returns structurei, rotation, basis, coeff instead of just structure.
% For compliance with the NRSFM full-initialization interface.

function [structure, rotations, basis, coeff] = ...
    nrsfm_init_find_structure_adaptor(f, projections, rotations, K)
  structure = f(projections, rotations);
  [basis, coeff] = factorize_structure(structure, K);
end
