function R = rotx(theta)
  c = cos(theta);
  s = sin(theta);
  R = [1, 0, 0; 0, c, -s; 0, s, c];
end
