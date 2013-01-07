COPYRIGHT © 2011-2012 Yuchao Dai, Hongdong Li, Mingyi He

   College of Engineering and Computer Science
   The Australian National University, Australia.
   School of Electronics and Information
   Northwestern Polytechnical University, P.R.China.

This code is an implementation of the methods described in:

      Yuchao Dai, Hongdong Li, Mingyi He: A simple prior-free method
      for non-rigid structure-from-motion factorization. CVPR 2012: 2018-2025

This software is distributed WITHOUT ANY WARRANTY. Use of this software is granted for research conducted at research institutions only. 
Commercial use of this software is not allowed. Corporations interested in the use of this software should contact the authors. 
If you use this code for a scientific publication, please reference the paper above.

USAGE:

Please see the demo files "NRSFM_example.m" for usage information. These scripts were tested with MATLAB versions R2010a.
Note that Yalmip and SDPT3 are used in the programs, which are available at http://users.isy.liu.se/johanl/yalmip/
and http://www.math.nus.edu.sg/~mattohkc/sdpt3.html

FEEDBACK:

Your feedback is greatly welcome. Please send bug reports, suggestions, and/or new results to:

   daiyuchao@gmail.com, hongdong.li@anu.edu.au

In the future, updated versions of this software will be available at:

   http://users.cecs.anu.edu.au/~yuchao/


The files in the directory includes:
ReadMe.txt:                                This file.
NRSFM_example.m:                           Example showing how to use this program.
NRSFM_BMM.m:                               Implementation of the methods.
Estimate_G_k.m:                            Estimate the corrective matrix G_k
Recover_Roation.m                          Recover the rotation
S_to_Shat.m                                Project the shape to K-shape bases
shape_recovery_fpca_s_sharp.m              Recover shape given measurements W and rotation
evalQ_regularization.m                     Evaluate the estimation of Q
Initialization.m                           Path setting with Yalmip and SDPT3, should be adapted to your setting
Rotate_Struct_Shape.m                      Output the rotated shape, modified from compareStructs.m

compareStructs.m                           These five files are part of the nrsfm_matlabCode downloaded from 
rotateStruct.m                             http://cvlab.lums.edu.pk/nrsfm/
procrist.m
findRotation.m
compareRotation.m           