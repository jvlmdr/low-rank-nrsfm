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
%function Shat = S_to_Shat(S,K)
%Transform S to make it low rank (3K) but each triple rows is of rank K
%Inputs:
%      S: 3F x P shape matrix
%      K: Number of shape basis
%Output:
%     Shat: 3F x P shape matrix satisfy the K-basis linear combination
%     constraint
%Reference: Yuchao Dai, Hongdong Li, Mingyi He: A simple prior-free method
%for non-rigid structure-from-motion factorization. CVPR 2012: 2018-2025
%Author: Yuchao Dai
%Contact Information: daiyuchao@gmail.com, yuchao.dai@anu.edu.au,
%hongdong.li@anu.edu.au
%Last update: 2012-10-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Shat = S_to_Shat(S,K)

[nPose nPts] = size(S);
nPose = nPose/3;

S_X = S(1:3:end,:);
S_Y = S(2:3:end,:);
S_Z = S(3:3:end,:);

S_Star = [S_X S_Y S_Z];

[U_S D_S V_S] = svd(S_Star);
D_S(K+1:end,K+1:end) = 0;
S_sharp = U_S*D_S*V_S';
 
Shat = zeros(3*nPose,nPts);

Shat(1:3:end,:) = S_sharp(:,1:nPts);
Shat(2:3:end,:) = S_sharp(:,nPts+1:2*nPts);
Shat(3:3:end,:) = S_sharp(:,2*nPts+1:3*nPts);
