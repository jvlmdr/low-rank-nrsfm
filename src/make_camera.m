function camera = make_camera(P)
  if nargin == 0
    P = [];
  end

  camera = struct('P', P);
end
