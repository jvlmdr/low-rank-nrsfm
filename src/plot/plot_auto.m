function plot_auto(ax, x, varargin)
  if length(size(x)) == 2
    if size(x, 2) == 2
      plot(ax, x(:, 1), x(:, 2), varargin{:});
    else
      plot3(ax, x(:, 1), x(:, 2), x(:, 3), varargin{:});
    end
  else
    if size(x, 2) == 2
      plot(ax, x(:, :, 1), x(:, :, 2), varargin{:});
    else
      plot3(ax, x(:, :, 1), x(:, :, 2), x(:, :, 3), varargin{:});
    end
  end
end
