%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code from http://cvlab.lums.edu.pk/nrsfm/
% Ijaz Akhter, Yaser Sheikh, Sohaib Khan, Takeo Kanade,
%  "Nonrigid Structure from Motion in Trajectory Space",
%  in Proceedings of 22nd Annual Conference on Neural Information Processing Systems, NIPS 2008, Vancouver, Canada, pp. 41-48, December 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Shat] = rotateStruct(Sh, Rs)

F = size(Rs,1)/2;

for i=1:F
    Y1 = Rs(2*i-1:2*i, :);
    Y1(3, :) = cross(Y1(1,:), Y1(2,:));
    Y2 = [Y1(1:2, :); -Y1(3,:)];
    if det(Y1)>0
        Y = Y1;
    else
        Y = Y2;
    end;
    Shat(3*i-2:3*i, :) = Y*Sh(3*i-2:3*i, :);
end;
