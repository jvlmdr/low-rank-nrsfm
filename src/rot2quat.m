function q = rot2quat(R)
  Qxx = R(1, 1);
  Qxy = R(1, 2);
  Qxz = R(1, 3);
  Qyx = R(2, 1);
  Qyy = R(2, 2);
  Qyz = R(2, 3);
  Qzx = R(3, 1);
  Qzy = R(3, 2);
  Qzz = R(3, 3);

  t = trace(R);
  r = sqrt(1 + t);
  w = 0.5 * r;
  x = sign(Qzy - Qyz) * 0.5 * sqrt(1 + Qxx - Qyy - Qzz);
  y = sign(Qxz - Qzx) * 0.5 * sqrt(1 - Qxx + Qyy - Qzz);
  z = sign(Qyx - Qxy) * 0.5 * sqrt(1 - Qxx - Qyy + Qzz);

  q = [w; x; y; z];

%  K = [Qxx - Qyy - Qzz,       Qyx + Qxy,       Qzx + Qxz,       Qyz - Qzy;
%             Qyx + Qxy, Qyy - Qxx - Qzz,       Qzy + Qyz,       Qzx - Qxz;
%             Qzx + Qxz,       Qzy + Qyz, Qzz - Qxx - Qyy,       Qxy - Qyx;
%             Qyz - Qzy,       Qzx - Qxz,       Qxy - Qyx, Qxx + Qyy + Qzz] / 3;
%
%  [V, D] = eig(K);
%  d = diag(D);
%  [~, ind] = sort(abs(d), 'descend');
%  d = d(ind)
%  V = V(:, ind);
%  q = V(:, 1);
end
