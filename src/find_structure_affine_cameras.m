% Finds S that minimizes
%   || reshape(S) ||_{*}
% subject to
%   w_ti = R_t s_ti.
%
% Parameters:
% W -- 2F x P
% R -- 2 x 3 x F

function X = find_structure_affine_cameras(W, Q, use_3P, settings)
  N = size(W, 2);
  F = size(Q, 3);

  % [2F, N] -> [2, F, N] -> [2, N, F]
  W = reshape(W, [2, F, N]);
  W = permute(W, [1, 3, 2]);

  converged = false;
  num_iter = 0;
  rho = settings.rho;

  X = zeros(3, N, F);
  Z = zeros(3, N, F);
  U = X - Z;

  while ~converged && num_iter < settings.max_iter
    prev_Z = Z;

    % Constrained least-squares.
    % Independent per projection, same matrix for each frame.
    V = Z - U;

    for t = 1:F
      P = eye(3);
      q = -V(:, :, t);
      A = Q(:, :, t);
      b = W(:, :, t);
      M = [P, A'; A, zeros(2, 2)];
      x = M \ [-q; b];
      X(:, :, t) = x(1:3, :);
    end

    % Singular value soft thresholding.
    V = X + U;
    if use_3P
      % [3, N, F] -> [3N, F]
      V = reshape(V, [3 * N, F]);
    else
      % [3, N, F] -> [3, F, N] -> [3F, N]
      V = reshape(permute(V, [1, 3, 2]), [3 * F, N]);
    end
    Z = singular_value_soft_threshold(V, 1 / rho);
    if use_3P
      % [3N, F] -> [3, N, F]
      Z = reshape(Z, [3, N, F]);
    else
      % [3F, N] -> [3, F, N] -> [3, N, F]
      Z = permute(reshape(Z, [3, F, N]), [1, 3, 2]);
    end

    % Update multipliers.
    R = X - Z;
    U = U + R;
    S = rho * (Z - prev_Z);

    norm_r = norm(R(:));
    norm_s = norm(S(:));

    fprintf('%12d %12g %12g %12g\n', num_iter, rho, norm_r, norm_s);

    if num_iter > 0 && num_iter < settings.max_iter / 2
      if norm_r ~= 0 && norm_s ~= 0
        if norm_r > settings.mu * norm_s
          rho = rho * settings.tau_incr;
        elseif norm_s > settings.mu * norm_r
          rho = rho / settings.tau_decr;
        end
      end
    end

    num_iter = num_iter + 1;
  end

  % [3, N, F] -> [3, F, N] -> [3F, N]
  X = reshape(permute(Z, [1, 3, 2]), [3 * F, N]);
end
