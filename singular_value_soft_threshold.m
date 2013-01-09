% Finds X that minimizes
%   1/2 || X - Y ||_F^2 + tau || X ||_*

function X = singular_value_soft_threshold(Y, tau)
  [m, n] = size(Y);
  %disp(['Computing ' num2str(m) ' x ' num2str(n) ' SVD...']);

  [U, S, V] = svd(Y, 'econ');
  s = diag(S);
  s = max(0, s - tau);
  X = U * diag(s) * V';
end
