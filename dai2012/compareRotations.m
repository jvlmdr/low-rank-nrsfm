%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code from http://cvlab.lums.edu.pk/nrsfm/
% Ijaz Akhter, Yaser Sheikh, Sohaib Khan, Takeo Kanade,
%  "Nonrigid Structure from Motion in Trajectory Space",
%  in Proceedings of 22nd Annual Conference on Neural Information Processing Systems, NIPS 2008, Vancouver, Canada, pp. 41-48, December 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errR, x, Rs1] = compareRotations(Rs, Rs1)

Rs1 = procrust(Rs, Rs1);
F = size(Rs, 1)/2;
for i=1:F
   errR(i) = norm(Rs1(2*i-1:2*i, :) - Rs(2*i-1:2*i, :));
end