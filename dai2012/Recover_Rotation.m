%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright © 2011,2012  Yuchao Dai, Hongdong Li, Mingyi He
%     This file is part of NRSFM_DLH.
% 
%     NRSFM_DLH is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     NRSFM_DLH is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with NRSFM_DLH.  If not, see <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [R_Recover Rsh] = Recover_Rotation(PI_Hat,G)
%Inputs:
%      PI_Hat: 2F*3K from the factorization of W
%      G: 3K*3 Corrective matrix
%Outputs:
%      R_Recover: Rotation matrix 2F*3
%      Rsh: Rotation matrix 2F*3F
%Reference: Yuchao Dai, Hongdong Li, Mingyi He: A simple prior-free method
%for non-rigid structure-from-motion factorization. CVPR 2012: 2018-2025
%Author: Yuchao Dai
%Contact Information: daiyuchao@gmail.com, yuchao.dai@anu.edu.au,
%hongdong.li@anu.edu.au
%Last update: 2012-10-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R_Recover Rsh] = Recover_Rotation(PI_Hat,G)

nPose = size(PI_Hat,1)/2;  %Number of frames

R_Recover = zeros(2*nPose,3); %Recovered matrix 2F*3

c_fk1 = zeros(nPose,1);
c_fk2 = zeros(nPose,1);

for f=1:nPose
    Eq = PI_Hat(2*f-1:2*f,:)*G*G'*PI_Hat(2*f-1:2*f,:)';

    c_fk1(f) = sqrt(abs(Eq(1,1))); %be determined by the rotation
    c_fk2(f) = sqrt(abs(Eq(2,2))); %be determined by the rotation
    
    R_Recover(2*f-1,:) = PI_Hat(2*f-1,:)*G/(c_fk1(f));
    R_Recover(2*f-0,:) = PI_Hat(2*f-0,:)*G/(c_fk2(f));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Here we have obtained the rotation for each frame and the coefficient,
%however there is a sign ambiguity, then the problem is how to fix it, the
%constraint between two frames is that the angle is less than 90 degrees.
r1 = R_Recover(1,:);
r2 = R_Recover(2,:);
r3 = cross(r1,r2);
R_f_1 = [r1;r2;r3];

if(det(R_f_1) < 0)
    R_f_1 = - R_f_1;
end

R_Recover(1:2,:) = R_f_1(1:2,:);

theta1 = zeros(nPose,1);
theta2 = zeros(nPose,1);

for f = 2:nPose
    r1 = R_Recover(2*f-1,:);
    r2 = R_Recover(2*f-0,:);
    r3 = cross(r1,r2);
    R_f = [r1;r2;r3];
    
    if(det(R_f) < 0)
        R_f = - R_f;
    end
        
    theta1(f) = acos((trace(R_f'*R_f_1)-1)/2)*180/pi;
    
    if(theta1(f) > 90)
        R_f(1:2,:) = - R_f(1:2,:);
    end
    
    theta2(f) = acos((trace(R_f'*R_f_1)-1)/2)*180/pi;
        
    R_Recover(2*f-1:2*f,:) = R_f(1:2,:);
    
    R_f_1 = R_f;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rsh is the big R 2F*3F diagonal
for i = 1:nPose
    Rsh(2*(i-1)+1:2*(i-1)+2,3*(i-1)+1:3*(i-1)+3) = R_Recover(2*i-1:2*i,:);
end