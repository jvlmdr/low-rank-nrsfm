% Finds S that minimizes
%   || reshape(S) ||_{*}
% subject to
%   w_i = R_i s_i    i = 1, ..., N.

function X = find_structure(projections, use_3P, settings)

  N = length(projections.tracks);
  F = projections.num_frames;

  converged = false;
  num_iter = 0;
  rho = settings.rho;

  X = zeros(3, F, N);
  Z = zeros(3, F, N);
  U = X - Z;

  while ~converged && num_iter < settings.max_iter
    prev_Z = Z;

    % Constrained least-squares, independent per projection.
    V = Z - U;
    % If there is no equation for x_{ti}, then set to v_{ti}.
    X = V;

    for i = 1:N
      track = projections.tracks(i);
      m = length(track.frames);

      for j = 1:m
        t = track.frames(j);
        P = eye(3);
        q = -V(:, t, i);
        A = track.equations(j).A;
        b = track.equations(j).b;
        M = [P, A'; A, zeros(2, 2)];
        x = M \ [-q; b];
        X(:, t, i) = x(1:3);
      end
    end

    % Singular value soft thresholding.
    V = X + U;
    if use_3P
      % Reshape to [3N, F] via [3, N, F].
      V = reshape(permute(V, [1, 3, 2]), [3 * N, F]);
    else
      % Reshape to [3F, N].
      V = reshape(V, [3 * F, N]);
    end
    Z = singular_value_soft_threshold(V, 1 / rho);
    if use_3P
      % Return to [3, F, N] via [3, N, F].
      Z = permute(reshape(Z, [3, N, F]), [1, 3, 2]);
    else
      % Return to [3, F, N].
      Z = reshape(Z, [3, F, N]);
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

  X = reshape(Z, [3 * F, N]);
end
