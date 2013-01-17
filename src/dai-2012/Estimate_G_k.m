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
%function G = Estimate_G_k(PI_Hat,V_A,count,K)
%Recover the corrective matrix G_k from trace norm minimization (Solved in
%SDP) and nonlinear optimization
%Inputs:
%      PI_Hat: 2F*3K  from factorization of measurement matrix W
%      V_A:    (3*K)*(3*K+1)/2  x 2K^2 - K.
%      K: Number of shape bases
%Output:
%      G:
%Reference: Yuchao Dai, Hongdong Li, Mingyi He: A simple prior-free method
%for non-rigid structure-from-motion factorization. CVPR 2012: 2018-2025
%Author: Yuchao Dai
%Contact Information: daiyuchao@gmail.com, yuchao.dai@anu.edu.au,
%hongdong.li@anu.edu.au
%Last update: 2012-10-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = Estimate_G_k(PI_Hat,V_A,count,K)

%Using nuclear norm to find the coefficients for the basis
X_Line = V_A(:,count-(2*K^2-K)+1:end); %This is the null space
    
t = sdpvar(2*K^2-K,1,'full');

Q_Line = X_Line*t;  %The constructed measurement matrix
    
%Transform the data to matrix from vector
Q_Hat = sdpvar(3*K,3*K,'full');
count = 0;
for i = 1:3*K
    for j = i:3*K
        count = count + 1;
        Q_Hat(i,j) = Q_Line(count);
        Q_Hat(j,i) = Q_Line(count);
    end
end

%Minimize Nuclear norm of W and Nuclear norm of lambda
XX = sdpvar(size(Q_Hat,1),size(Q_Hat,1));
YY = sdpvar(size(Q_Hat,2),size(Q_Hat,2));  % auxiliary matrices to solve the sdp

Z = [XX Q_Hat
     Q_Hat' YY];      % the constructed matrix

b = ones(2*K^2-K,1);

Fun = [Z>=0, b'*t==1]; %  t(1) == 1 norm(t)==1  ,Q_Hat1>=0 t(1) == -1  b'*t==1
h = 1/2*(trace(XX) + trace(YY));

solvesdp(Fun,h);
    
scale = b'*t;

Q_Hat = double(Q_Hat/scale);

%--------------------------------------------------------------------------
% for k = 1:100
%     [U_Hat D_Hat V_Hat] = svd(Q_Hat);
% D_Hat(3,3)/D_Hat(4,4)
% D_Hat(4:end,4:end) = 0;
% Q_Hat = U_Hat*D_Hat*V_Hat;  %Rank-3
% 
% count = 0;
% for i = 1:3*K
%     for j = i:3*K
%         count = count + 1;
%         q_Line(count,1) = Q_Hat(i,j);
%     end
% end
% 
% %Projection on subspace
% coefficient = pinv(X_Line)*q_Line;
% 
% q_Line = X_Line*coefficient/norm(coefficient);
% 
% count = 0;
% for i = 1:3*K
%     for j = i:3*K
%         count = count + 1;
%         Q_Hat(i,j) = q_Line(count);
%         Q_Hat(j,i) = q_Line(count);
%     end
% end
% end
%--------------------------------------------------------------------------


[U_Hat D_Hat V_Hat] = svd(Q_Hat);

G_Hat = U_Hat(:,1:3)*sqrt(D_Hat(1:3,1:3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nonlinear optimization
disp('  Nonlinear refinement of G...(wait).'); 

options = optimset('Display', 'Final', 'Diagnostics','off','Largescale', 'off', 'MaxFunEval',200000,'MaxIter',5000,'TolFun',1e-10,'TolX',1e-10);

 [G, fval] = fminunc(@evalQ_regularization,G_Hat,options,PI_Hat);  %fminunc  fminsearch

%[G, fval] = fminsearch(@evalQ_regularization,G_Hat,options,PI_Hat);
