function plot_camera(ax, P, axis_length)
  R = P(:, 1:3);
  d = P(:, 4);
  c = -R' * d;

  line_width = 1;

  i = [c, c + axis_length * R(1, :)'];
  j = [c, c + axis_length * R(2, :)'];
  k = [c, c + axis_length * R(3, :)'];

  plot3(i(1, :), i(3, :), i(2, :), 'r-', 'LineWidth', line_width);
  plot3(j(1, :), j(3, :), j(2, :), 'g-', 'LineWidth', line_width);
  plot3(k(1, :), k(3, :), k(2, :), 'b-', 'LineWidth', line_width);
end
