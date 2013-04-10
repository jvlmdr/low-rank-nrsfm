function Q = find_corrective_gram_matrix_sweep(A, c)
  lambdas = logspace(-6, 12, 19)';
  num_lambdas = length(lambdas);

  m = size(A, 1);
  n = size(A, 2);
  % p (p + 1) = 2n
  % p < sqrt(2n) < p + 1
  p = floor(sqrt(2 * n));
  [~, subset] = construct_symmetric(p);

  min_e = inf;
  Q = [];
  arg = [];

  for i = 1:num_lambdas
    Q_i = find_corrective_gram_matrix(A, c, lambdas(i));

    Q_low_rank = project_psd_rank(Q_i, 3);
    q = Q_low_rank(subset)';

    e = norm(A * q);
    if e < min_e
      Q = Q_i;
      arg = i;
      min_e = e;
    end

    fprintf('lambda: % .2e  e: % .2e\n', lambdas(i), e);
  end

  fprintf('Best solution used lambda = %g\n', lambdas(arg));
end
