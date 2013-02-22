% Parameters:
% S, X -- 3F x P non-rigid shape

function residual = min_shape_error(S, X)
  F = size(S, 1) / 3;
  P = size(S, 2);

  % [3F, P] -> [3, F, P] -> [P, 3, F]
  S = reshape(S, [3, F, P]);
  S = permute(S, [3, 1, 2]);
  X = reshape(X, [3, F, P]);
  X = permute(X, [3, 1, 2]);

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
