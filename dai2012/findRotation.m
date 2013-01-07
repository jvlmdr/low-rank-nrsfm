%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code from http://cvlab.lums.edu.pk/nrsfm/
% Ijaz Akhter, Yaser Sheikh, Sohaib Khan, Takeo Kanade,
%  "Nonrigid Structure from Motion in Trajectory Space",
%  in Proceedings of 22nd Annual Conference on Neural Information Processing Systems, NIPS 2008, Vancouver, Canada, pp. 41-48, December 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y] = findRotation(S, Sh)

F = size(S, 1)/3;
S1 = [];
S2 = [];
for i=1: F
    S1 = [S1, S(3*i-2:3*i, :)];
    S2 = [S2, Sh(3*i-2:3*i, :)];
end;
Y = S1/S2;
% Y = inv(S2*S2')*S1*S2';
[U, D, V] = svd(Y);
Y = U*V';