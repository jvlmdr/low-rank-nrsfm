% Find G = Q Q', where Q is the 3K x 3 matrix which best corrects P_hat.
% Not exactly Dai's method since it uses ADMM instead of SDP solver.

function R = find_rotations_dai(P_hat)

%  F = size(P_hat, 1) / 2;
%  K = size(P_hat, 2) / 3;
%
%  H = construct_symmetric(3 * K);
%  n = (3 * K) * (3 * K + 1) / 2;
%
%  % Get constraints that make A orthogonal.
%  [A, c] = rotation_constraints(P_hat);
%  m = size(A, 2);
%  % c' x = d.
%  C = c';
%  d = 1;
%
%  % Find an estimate using trace-norm minimization.
%  G = find_corrective_gram_matrix(A, c, lambda);
%  g = H \ G(:);
%  fprintf('rank(G) = %d\n', rank(G));
%  fprintf('|| A g || / || g || = %g\n', norm(A * g) / norm(g));
%
%  [U D V] = svd(G);
%  Q_hat = U(:,1:3)*sqrt(D(1:3,1:3));
%  disp('  Nonlinear refinement of G...(wait).'); 
%  options = optimset('Display', 'Final', 'Diagnostics','off','Largescale', 'off', 'MaxFunEval',200000,'MaxIter',5000,'TolFun',1e-10,'TolX',1e-10);
%  [Q, fval] = fminunc(@evalQ_regularization,Q_hat,options,P_hat);

  Q = find_corrective_matrix_dai(P_hat);
  R = Recover_Rotation(P_hat, Q);
end
