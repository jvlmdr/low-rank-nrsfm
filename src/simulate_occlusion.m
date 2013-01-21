function mask = simulate_occlusion(num_frames, num_points, missing, ...
    mean_length)
  p_occlude = 1 / mean_length;
  p_unocclude = (1 - missing) / missing * p_occlude;

  % Probability of transitioning visible = 0 to 1 and 1 to 0, respectively.
  p_change = [p_unocclude; p_occlude];
  % Initial state.
  p_visible = 1 - missing;
  mask = zeros(num_frames, num_points);

  for i = 1:num_points
    r_init = rand();
    r = rand(num_frames - 1, 1);

    % Frames t:num_frames still need to be labelled.
    t = 1;
    v = r_init < p_visible;
    visible = zeros(0, 1);

    while t <= num_frames
      % Look for the next transition event from state v.
      j = find(r(t:end) < p_change(v + 1), 1);

      if isempty(j)
        % No changes until the end.
        j = num_frames - t + 1;
      end

      % Next j frames are the same.
      visible(t:t + j - 1) = v;
      t = t + j;
      v = 1 - v;
    end

    mask(:, i) = visible;
  end
end
