function [mask, index] = separate_occlusions(original)
  num_frames = size(original, 1);
  num_points = size(original, 2);

  mask = sparse(num_frames, 0);

  index = [];
  num_tracks = 0;
  for i = 1:num_points
    % Still need to find visible segments in frames t..num_frames.
    t = 1;
    found_all = false;

    while ~found_all
      % Find the next visible frame.
      j = find(original(t:num_frames, i), 1);

      if isempty(j)
        found_all = true;
      else
        t = t + j - 1;

        % There is a visible segment.
        % Find the next occluded frame.
        j = find(~original(t:num_frames, i), 1);

        if isempty(j)
          % Go to end of sequence.
          j = num_frames - t + 1;
        end

        % Extract.
        visible = sparse(num_frames, 1);
        visible(t:t + j - 1) = 1;
        mask = [mask, visible];
        index = [index, i];
        num_tracks = num_tracks + 1;
        t = t + j;
      end
    end
  end
end
