% G = refine_corrective_triple_nonlinear(M_hat, G_init)
%
% Parameters:
% M_hat -- 2 x 3K x F
% G_init -- 3K x 3
%
% Returns:
% G -- 3K x 3

function G = refine_corrective_triple_nonlinear(M_hat, G, max_iter, tol)
  assert(ndims(M_hat) == 3);
  assert(size(M_hat, 1) == 2);
  K = size(M_hat, 2) / 3;
  F = size(M_hat, 3);

  assert(ndims(G) == 2);
  assert(size(G, 1) == 3 * K);
  assert(size(G, 2) == 3);

  [A, c] = rotation_constraints(projections_to_matrix(M_hat));
  [H, subset] = construct_symmetric(3 * K);

%  % Rescale current solution.
%  G = G / norm(G(:));
%  Q = G * G';
%  q = Q(subset)';
%  r = 1/2 * norm(A * q) ^ 2;
%  fprintf('Initial residual: %g\n', r);
%  r_init = r;

  G = 1 / norm(G(:)) * G;

  % Solve.
  if exist('find_structure_low_rank_nonlinear_mex', 'file')
    G = refine_corrective_triple_nonlinear_mex(M_hat, G, max_iter, tol);
  else
    error('Not compiled');
  end

%  G = G / norm(G(:));
%  Q = G * G';
%  q = Q(subset)';
%  r = 1/2 * norm(A * q) ^ 2;
%  fprintf('Initial residual: %g\n', r_init);
%  fprintf('Final residual: %g\n', r);
end
