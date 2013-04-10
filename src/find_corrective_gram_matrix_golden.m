function Q = find_corrective_gram_matrix_golden(F, E, tol)
  a = log(1e-6);
  b = log(1e9);

  m = size(F, 1);
  n = size(F, 2);
  % p (p + 1) = 2n
  % p < sqrt(2n) < p + 1
  p = floor(sqrt(2 * n));
  [~, subset] = construct_symmetric(p);

  phi = (sqrt(5) - 1) / 2;

  min_e = inf;
  Q = [];
  arg = [];

  Q = find_corrective_gram_matrix(F, E, exp(a));
  Q = project_psd_rank(Q, 3);
  e_a = norm(F * Q(subset)');

  Q = find_corrective_gram_matrix(F, E, exp(b));
  Q = project_psd_rank(Q, 3);
  e_b = norm(F * Q(subset)');

  c = a + phi * (b - a);
  Q = find_corrective_gram_matrix(F, E, exp(c));
  Q_c = Q;
  Q = project_psd_rank(Q, 3);
  e_c = norm(F * Q(subset)');

  converged = false;

  while ~converged
    % Let d be the new parameter value.
    if c < (a + b) / 2
      d = a + phi * (b - a);
    else
      d = a + (1 - phi) * (b - a);
    end

    Q = find_corrective_gram_matrix(F, E, exp(d));
    Q_d = Q;
    Q = project_psd_rank(Q, 3);
    e_d = norm(F * Q(subset)');

    % Let x be the left value, y be the right value.
    x = c;
    e_x = e_c;
    y = d;
    e_y = e_d;
    if y < x
      tmp = x; x = y; y = tmp;
      tmp = e_x; e_x = e_y; e_y = tmp;
    end

    % Let c be the min value.
    if e_d < e_c
      c = d;
      e_c = e_d;
      Q_c = Q_d;
    end

    fprintf(' a: %8.4f   x: %8.4f   y: %8.4f   b: %8.4f\n', a, x, y, b);
    fprintf('fa: %8.4f  fx: %8.4f  fy: %8.4f  fb: %8.4f\n', e_a, e_x, e_y, e_b);

    if e_x > e_a || e_y > e_b
      warning('Function is not unimodal');
    end

    if e_x < e_y
      % Keep a, move b in.
      b = y;
      e_b = e_y;
    else
      % Keep b, move a in.
      a = x;
      e_a = e_x;
    end

    if abs(x - y) <= 1e-6 + tol * max(abs(x), abs(y))
      converged = true;
    end
  end

  Q = Q_c;
  fprintf('lambda: %g\n', exp(c));
end
