function R = quat2rot(q)
  w = q(1);
  x = q(2);
  y = q(3);
  z = q(4);

  ww = w * w;
  wx = w * x;
  wy = w * y;
  wz = w * z;
  xx = x * x;
  xy = x * y;
  xz = x * z;
  yy = y * y;
  yz = y * z;
  zz = z * z;

  R = [1 - 2 * yy - 2 * zz,     2 * xy - 2 * wz,     2 * xz + 2 * wy;
           2 * xy + 2 * wz, 1 - 2 * xx - 2 * zz,     2 * yz - 2 * wx;
           2 * xz - 2 * wy,     2 * yz + 2 * wx, 1 - 2 * xx - 2 * yy];

%  a = q(1);
%  b = q(2);
%  c = q(3);
%  d = q(4);
%
%  aa = a * a;
%  ab = a * b;
%  ac = a * c;
%  ad = a * d;
%  bb = b * b;
%  bc = b * c;
%  bd = b * d;
%  cc = c * c;
%  cd = c * d;
%  dd = d * d;
%
%  R(1, 1) = aa + bb - cc - dd;
%  R(1, 2) = 2 * (bc - ad);
%  R(1, 3) = 2 * (ac + bd);
%  R(2, 1) = 2 * (ad + bc);
%  R(2, 2) = aa - bb + cc - dd;
%  R(2, 3) = 2 * (cd - ab);
%  R(3, 1) = 2 * (bd - ac);
%  R(3, 2) = 2 * (ab + cd);
%  R(3, 3) = aa - bb - cc + dd;

  R = R / norm(q)^2;
end
