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
%[f] = evalQ_regularization(q,PI_Hat)
%Evaluate the estimation of q
%Input:
%     q: Current estimation of q with size 3K*3
%     PI_Hat:Measurements with size 2F*3K
%Output:
%     f: Characterizing the difference from orthogonal constraints
%Reference: Yuchao Dai, Hongdong Li, Mingyi He: A simple prior-free method
%for non-rigid structure-from-motion factorization. CVPR 2012: 2018-2025
%Author: Yuchao Dai
%Contact Information: daiyuchao@gmail.com, yuchao.dai@anu.edu.au,
%hongdong.li@anu.edu.au
%Last update: 2012-10-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = evalQ_regularization(q,PI_Hat)

Rbar = PI_Hat*q;        % Rbar(2F*3), rotation matrices

Rbar1 = Rbar(1:2:end, :);
Rbar2 = Rbar(2:2:end, :);

F = size(Rbar,1)/2;
C = zeros(1,F);
diffArr = zeros(1,F);

zerosArr = sum(Rbar1.*Rbar2,2);

for i = 1:F
    C(i) = sum(Rbar1(i,:).^2);
    diffArr(i) = 1 - sum(Rbar2(i,:).^2)/C(i);
    
    zerosArr(i) = zerosArr(i)/C(i);
end

f = sum(diffArr.^2) + 4*sum(zerosArr.^2);