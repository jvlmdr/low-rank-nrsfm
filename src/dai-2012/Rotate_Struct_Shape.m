%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code Modfied from the function  compareStructs
% from http://cvlab.lums.edu.pk/nrsfm/
% Ijaz Akhter, Yaser Sheikh, Sohaib Khan, Takeo Kanade,
%  "Nonrigid Structure from Motion in Trajectory Space",
%  in Proceedings of 22nd Annual Conference on Neural Information Processing Systems, NIPS 2008, Vancouver, Canada, pp. 41-48, December 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Shat = Rotate_Struct_Shape(S, Shat)

s = mean(std(S, 1, 2));
S = S/s;
sm = mean(S,2);
S = S - sm*ones(1,size(S,2));

Shat = Shat/s;
F = size(S,1)/3;
P = size(S,2);

Y = findRotation(S, Shat);         % Procrust Alignment
for i=1:F
    Shat(3*i-2:3*i,:) = Y*Shat(3*i-2:3*i,:);
end