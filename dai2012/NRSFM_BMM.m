%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright ?2011,2012  Yuchao Dai, Hongdong Li, Mingyi He
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
%[Shat_BMM Shat_PI Shat_Smooth Rsh R_Recover Shape_Err_BMM Shape_Err_PI
%Shape_Err_Smooth Shape_Err_SVP Shape_BMM_Time err3D Rotation_Err] = NRSFM_BMM(W,K,rotStruct,S_GT,Rs)
%Algorithm implementation of our rotation and shape recovery, including SDP
%based rotation recovery+nonlinear optimization+pseudo-inverse method +
%block matrix method+ smoothness constrained method
%Inputs:
%        W: Measurement matrix with dimension 2F * P
%        K: Estimation of number of shape bases
%        rotStruct: Boolean value
%        S_GT: Ground truth deformable shape 3F * P
%        Rs: Ground truth rotation           2F * 3
%Outputs:
%        Shat_BMM: Deformable shape from block matrix method
%        Shat_PI:  Deformable shape from pseudo-inverse method
%        Shat_Smooth: Deformable shape from smoothness constrained method
%        Rsh: Recovered rotation matrix 2F * 3F
%        R_Recover: Recovered rotation matrix 2F * 3
%        Shape_Err_BMM: Shape recovery error of block matrix method
%        Shape_Err_PI:  Shape recovery error of pseudo-inverse method
%        Shape_Err_Smooth: Shape recovery error of smoothness constrained
%        method
%        Rotation_Err: Rotation recovery error
%Reference: Yuchao Dai, Hongdong Li, Mingyi He: A simple prior-free method
%for non-rigid structure-from-motion factorization. CVPR 2012: 2018-2025
%Author: Yuchao Dai
%Contact Information: daiyuchao@gmail.com, yuchao.dai@anu.edu.au,
%hongdong.li@anu.edu.au
%Last update: 2012-10-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Shat_BMM Shat_PI Shat_Improved Rsh R_Recover Shape_Err_BMM Shape_Err_PI Shape_Err_Smooth Rotation_Err] = NRSFM_BMM(W,K,rotStruct,S_GT,Rs)

nPose = size(W,1)/2;     %Number of frames
nPts = size(W,2);        %Number of points

% Normalize the S and W matrices when there is ground truth shape
if exist('S_GT','var')       
    s = mean(std(S_GT, 1, 2));
    S_GT = S_GT/s;
    sm = mean(S_GT,2);
    S_GT = S_GT - sm*ones(1,size(S_GT,2));
    W = W/s;
end

%--------------------------------------------------------------------------
%Centeralize the measurements
wm = mean(W,2);
W = W - wm*ones(1,size(W,2));

if (2*nPose>nPts)
    [U,D,V] = svd(W,0);
else
    [V,D,U] = svd(W',0);
end;

%--------------------------------------------------------------------------
%Fix the sign of U1
for i = 1:size(U,2)
    if(U(1,i) < 0)
        U(:,i) = -U(:,i);
    end
end

%--------------------------------------------------------------------------
%Obtain the initial factorization result
PI_Hat = U(:,1:3*K)*sqrt(D(1:3*K,1:3*K));

%Solve all the G_k or Q_k
A = zeros(2*nPose,(3*K)*(3*K));
A_Hat = zeros(2*nPose,(3*K)*(3*K+1)/2);

%Using the rotation constraint only    
%For each frame, there are two constraints
for f=1:nPose
    PI_Hat_f = PI_Hat(2*f-1:2*f,:);
    %Here we use the relationship that vec(ABA^T) = (A\odotA)vec(B)
    AA = kron(PI_Hat_f,PI_Hat_f);
        
    A(2*f-1,:) = AA(1,:) - AA(4,:);
    A(2*f-0,:) = AA(2,:);
        
    count = 0;
    for i=1:3*K
        for j=i:3*K
            count = count+1;
            if(i==j)
                A_Hat(2*f-1,count)=A(2*f-1,(3*K)*(i-1)+j);
                A_Hat(2*f-0,count)=A(2*f-0,(3*K)*(i-1)+j);
            else
                A_Hat(2*f-1,count)=A(2*f-1,(3*K)*(i-1)+j)+A(2*f-1,(3*K)*(j-1)+i);
                A_Hat(2*f-0,count)=A(2*f-0,(3*K)*(i-1)+j)+A(2*f-0,(3*K)*(j-1)+i);
            end
        end
    end
end

[U_A D_A V_A] = svd(A_Hat);

%--------------------------------------------------------------------------
%Fix the sign of V_A1: hongdong . 
for i = 1:size(V_A,2)
    if(V_A(1,i) < 0)
        V_A(:,i) = -V_A(:,i);
    end
end
%--------------------------------------------------------------------------

%Estimate G through SDP and nonlinear optimization
disp('Solve the first SDP to find G matrix... (wait)..'); 
G = Estimate_G_k(PI_Hat,V_A,count,K); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Recover the rotations, R_Recover 2F*3, Rsh 2F*3F

[R_Recover Rsh] = Recover_Rotation(PI_Hat,G);
%--------------------------------------------------------------------------
%If there exists ground truth rotation, evaluate the performance of
%rotation estimation, we use the same meteric as NIPS 20008
if exist('Rs','var')
    [errR] = compareRotations(Rs, R_Recover);
    disp('Rotation Estimation Error=')
    mean(errR),
    Rotation_Err = mean(errR);
else
    Rotation_Err = [];
end
%=========================================================================
%Recovery the deformable shape

%Method 1: Pseudo-inverse method
S_PI = pinv(Rsh)*W;                       % Rank-3K S
Shat_PI = S_to_Shat(S_PI,K);              % Transform S_PI to satisfy the K basis constraint

%Evaluate the shape recovery performance
if exist('S_GT','var')
    if rotStruct == 1
        Shat_PI = rotateStruct(Shat_PI, R_Recover);
    end

    errS = compareStructs(S_GT,Shat_PI);
      disp('Normalized 3D shape error by pseudoinverse method =');
    Shape_Err_PI = mean(errS), 
    
    %Shat_PI1 = Rotate_Struct_Shape(S_GT, Shat_PI1);
else
    Shape_Err_PI = NaN;
end

%--------------------------------------------------------------------------
%Method 2: Regularization constrained shape recovery
%S=(R^TR + \lambda H^TH)^+ R^T W with closed-form solution

%First order regular pattern
H = zeros(3*nPose);
for i = 1:3*nPose-3
    H(i,i) = -1;
    H(i,i+3) = 1;
end

%Trade-off parameter
lambda = 1e-3; %e-3;

%Closed-form solution, rank-3K as W
S_Improved = pinv(Rsh'*Rsh + lambda*H'*H)*Rsh'*W;

% Transform to K basis form
Shat_Improved = S_to_Shat(S_Improved,K);            

%Evaluate the performance
if exist('S_GT','var')
    if rotStruct == 1
        Shat_Improved = rotateStruct(Shat_Improved, R_Recover);
    end

    errS = compareStructs(S_GT,Shat_Improved);
    disp('Normalized 3D shape error by the improved pseudoinverse method =');%hongdong add this. 
    
    Shape_Err_Smooth = mean(errS),
    
    Shat_Improved = Rotate_Struct_Shape(S_GT, Shat_Improved);
else
    Shape_Err_Smooth = NaN;
end
%--------------------------------------------------------------------------
%Method 3: Block matrix method solved with fixed point continuation

disp('Start to solve shape matrix by the block matrix method using SDP/FPC ... (SDP/FPC is a rather slow process, be patient please....:-) hongdong ');
    

S_BMM = shape_recovery_fpca_s_sharp(W,Rsh,Shat_PI,K);

Shat_BMM = S_to_Shat(S_BMM,K);                   % Transform to K basis form

%--------------------------------------------------------------------------
%Evaluate the performance of BMM
if exist('S_GT','var')
    if rotStruct == 1
        [Shat_BMM] = rotateStruct(Shat_BMM, R_Recover);
    end

    errS = compareStructs(S_GT,Shat_BMM);
       disp('Normalized 3D shape error by block matrix method =');
    
    Shape_Err_BMM = mean(errS),
    
    Shat_BMM = Rotate_Struct_Shape(S_GT, Shat_BMM);
else
    Shape_Err_BMM = NaN;
end

if exist('S_GT','var')
    Shat_PI = Rotate_Struct_Shape(S_GT, Shat_PI);
end