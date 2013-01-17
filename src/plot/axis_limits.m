function limits = axis_limits(points)

lower = min(points, [], 1);
upper = max(points, [], 1);

limits = [lower; upper];
limits = limits(:)';

end
