% Parameters:
% S, X -- 3 x P x F non-rigid shape
% S is the reference, to which X is aligned.

function residual = min_shape_error(S, X)
  F = size(S, 1);
  P = size(S, 2);

  % [3, P, F] -> [P, 3, F]
  S = permute(S, [2, 1, 3]);
  X = permute(X, [2, 1, 3]);

  residual = 0;
  total = 0;

  for t = 1:F
    X_t = X(:, :, t);
    S_t = S(:, :, t);

    % Remove mean.
    X_t = bsxfun(@minus, X_t, mean(X_t));
    S_t = bsxfun(@minus, S_t, mean(S_t));

    % Normalize scale.
    x = sqrt(mean(sum(X_t .^ 2, 2), 1));
    s = sqrt(mean(sum(S_t .^ 2, 2), 1));
    X_t = s / x * X_t;

    % Align each frame individually.
    R = procrustes(X_t, S_t);
    D = S_t - X_t * R;
    residual = residual + mean(sqrt(sum(D .* D, 2)));
    total = total + mean(sqrt(sum(S_t .* S_t, 2)));
  end

  residual = residual / total;
end
