function [A, A_Hat, C] = rotation_constraints(PI_Hat)

% PI_Hat is 2F x 3K.
nPose = size(PI_Hat, 1) / 2;
K = size(PI_Hat, 2) / 3;

%Solve all the G_k or Q_k
A = zeros(2*nPose,(3*K)*(3*K));
A_Hat = zeros(2*nPose,(3*K)*(3*K+1)/2);

C = zeros(1, (3 * K) * (3 * K));

%Using the rotation constraint only    
%For each frame, there are two constraints
for f=1:nPose
    PI_Hat_f = PI_Hat(2*f-1:2*f,:);
    %Here we use the relationship that vec(ABA^T) = (A\odotA)vec(B)
    AA = kron(PI_Hat_f,PI_Hat_f);
        
    A(2*f-1,:) = AA(1,:) - AA(4,:);
    A(2*f-0,:) = AA(2,:);

    C = C + AA(1, :) + AA(4, :);
        
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

end
