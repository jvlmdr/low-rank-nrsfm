function trajectory = generate_velocity(num_frames, num_dims)
  trajectory = randn(num_frames - 1, num_dims);
  trajectory = [zeros(1, num_dims); cumsum(trajectory, 1)];
  trajectory = trajectory - ones(num_frames, 1) * mean(trajectory, 1);
end
