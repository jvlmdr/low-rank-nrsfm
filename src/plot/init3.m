function init(fig, lim)
  lim = lim([1, 2, 5, 6, 3, 4]);
  ax = gca(fig);
  axis(ax, lim);
  axis(ax, 'equal');
  axis(ax, 'manual');
  axis(ax, 'vis3d');
  hold(ax, 'on');
  grid(ax, 'on');
  box(ax, 'off');
  set(ax, 'YDir', 'reverse');

  xlabel(ax, 'x');
  ylabel(ax, 'z');
  zlabel(ax, 'y');

  set(ax, 'XTickLabel', {});
  set(ax, 'YTickLabel', {});
  set(ax, 'ZTickLabel', {});
end
