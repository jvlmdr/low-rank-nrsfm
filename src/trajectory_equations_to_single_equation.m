function [A, b] = trajectory_equations_to_single_equation(systems, num_frames)
  N = length(systems);
  F = num_frames;

  n = 3 * F * N;
  A = sparse(0, n);
  b = zeros(0, 1);

  for i = 1:N
    Q = systems(i).A;
    u = systems(i).b;
    m = size(Q, 1);
    A_i = sparse(m, 3 * F * N);
    A_i(:, (3 * F * (i - 1) + (1:3 * F))) = Q;
    A = [A; A_i];
    b = [b; u];
  end
end
