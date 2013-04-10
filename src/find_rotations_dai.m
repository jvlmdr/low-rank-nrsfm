% rotations = find_rotations_dai(projections, K)
%
% Parameters:
% projections -- 2 x P x F
% K -- Number of bases
%
% Returns:
% rotations -- 2 x 3 x F

function rotations = find_rotations_dai(projections, K)
  F = size(projections, 3);
  P = size(projections, 2);

  W = projections_to_matrix(projections);
  nPose = F;
  nPts = P;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Start of excerpt from NRSFM_BMM.m
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % End of excerpt from NRSFM_BMM.m
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  rotations = unstack_cameras(R_Recover);
end
