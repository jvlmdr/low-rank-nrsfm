function R = rotz(theta)
  c = cos(theta);
  s = sin(theta);
  R = [c, -s, 0; s, c, 0; 0, 0, 1];
end
