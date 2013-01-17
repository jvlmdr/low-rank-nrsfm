function init2(fig, lim)
  ax = gca(fig);
  axis(ax, lim);
  axis(ax, 'equal');
  axis(ax, 'manual');
  hold(ax, 'on');
  grid(ax, 'on');
  box(ax, 'off');
end
