% Returns structure and rotation instead of just structure.
% For compliance with the NRSFM interface.

function [structure, rotations] = nrsfm_find_structure_adaptor(f, ...
    projections, rotations)
  structure = f(projections, rotations);
end
