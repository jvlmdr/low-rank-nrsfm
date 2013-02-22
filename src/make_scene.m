function scene = make_scene(points, cameras, projections)
  if nargin == 0
    points = [];
    cameras = make_camera();
    projections = [];
  end

  scene = struct('points', points, 'cameras', cameras, ...
      'projections', projections);
end
