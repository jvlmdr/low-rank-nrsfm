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
%     Example use of the NRSFM method proposed in:  
%     Reference: Yuchao Dai, Hongdong Li, Mingyi He: A simple prior-free method
%     for non-rigid structure-from-motion factorization. CVPR 2012: 2018-2025
%     In this NRSFM_example, we load the image measurment matrix W, set
%     number of shape bases K. rotStruct is used to index whether we have
%     to rotate the shape and compare them. S gives the ground turth shape 
%     if available, Rs gives the ground truth rotation if available. 
%     In using your own datasets, you have to set the complete image
%     measurement matrix W, number of shape bases K and S, Rs, i.e. the
%     ground truth shape and rotation if available. 


%     All the datasets are tested. the results reported on CVPR12 paper were based on Matlab 7.10.0 (R2010a).
%     The following results are produced under Matlab R2011b and R2012a.  Note: the results are slightly different, which is to do with     
%     Matlab's SVD function's sign underterminacy. So far we still do nto have good method to fix it.  
% 
%--------------------------------------------------------------------------
%     Author: Yuchao Dai/hongdong Li 
%     Contact Information: daiyuchao@gmail.com, yuchao.dai@anu.edu.au,
%     hongdong.li@anu.edu.au
%     Last update: 2012-10-01
%     Note that the results may be slightly different on different machines and Matlab versions.
%     
%     All the datasets are downloaded from http://cvlab.lums.edu.pk/nrsfm/ and 
%     http://cbcsl.ece.ohio-state.edu/downloads.html

%     The hand-picked optimal K Values for different datasets are:

%     yoga            10     rotStruct = 0  Shape_Err_BMM 0.1150 Rotation_Err 0.0883
%     pickup          12     rotStruct = 0  Shape_Err_BMM 0.1731 Rotation_Err 0.1210
%     stretch         11     rotStruct = 0  Shape_Err_BMM 0.1034 Rotation_Err 0.0676
%     drink           12     rotStruct = 0  Shape_Err_BMM 0.0266 Rotation_Err 0.0071

%     dance           10     rotStruct = 1  Shape_Err_BMM 0.1864


%     face2           7      rotStruct = 1  Shape_Err_BMM 0.0303
%     walking2        8      rotStruct = 1  Shape_Err_BMM 0.1298
%     shark2          4      rotStruct = 1  Shape_Err_BMM 0.2357%--------------------------------------------------------------------------

clear
clc

load ../../data/akhter-2008/yoga

K = 10;

rotStruct = 0;

disp('===Nonrigid Factorization by Dai and Li====')
disp('=== Three algorithms  ( pseudo inverse, improved pseudo-inverse, and black-matrix-method) are provided.===');
disp('Note: the numerical results may be slightly different using different Matlab SVD function implementation.') ;


if exist('S','var')
    if exist('Rs','var')       % Normalize the S and W matrices
        [Shat_BMM Shat_PI Shat_Improve Rsh R_Recover Shape_Err_BMM Shape_Err_PI Shape_Err_Smooth Rotation_Err] = NRSFM_BMM(W,K,rotStruct,S,Rs);
    else
        [Shat_BMM Shat_PI Shat_Improve Rsh R_Recover Shape_Err_BMM Shape_Err_PI Shape_Err_Smooth Rotation_Err] = NRSFM_BMM(W,K,rotStruct,S);
    end
else
    [Shat_BMM Shat_PI Shat_Improve Rsh R_Recover Shape_Err_BMM Shape_Err_PI Shape_Err_Smooth Rotation_Err] = NRSFM_BMM(W,K,rotStruct);
end

