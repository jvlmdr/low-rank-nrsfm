% Parameters:
% S, X -- F x P x 3 non-rigid shape

function residual = min_shape_error(S, X)
  F = size(S, 1);
  P = size(S, 2);

  % [F, P, 3] -> [P, 3, F]
  S = permute(S, [2, 3, 1]);
  X = permute(X, [2, 3, 1]);

  residual = 0;
  total = 0;

  for t = 1:F
    % Align each frame individually.
    R = procrustes(X(:, :, t), S(:, :, t));
    D = S(:, :, t) - X(:, :, t) * R;
    residual = residual + mean(sqrt(sum(D .* D, 2)));
    total = total + mean(sqrt(sum(S(:, :, t) .* S(:, :, t), 2)));
  end

  residual = residual / total;
end
