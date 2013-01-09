function G = find_corrective_matrix_hongdong(PI_Hat)

% PI_Hat is 2F x 3K.
nPose = size(PI_Hat, 1) / 2;
K = size(PI_Hat, 2) / 3;

[~, A_Hat] = rotation_constraints(PI_Hat);
count = (3 * K) * (3 * K + 1) / 2;

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

end
