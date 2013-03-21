residuals = permute(residuals, [1, 3, 2]);

for i = 1:num_noises
  noise_stddev = noise_stddevs(i);
  figure;
  boxplot(residuals(:, :, i));
  hold on;
  plotSpread(residuals(:, :, i));
  hold off;
  grid on;
  axis([0, num_solvers + 1, 0, max(residuals(:)) * 1.05]);
  title(sprintf('Noise level %g', noise_stddev));
  set(gca, 'XTick', 1:num_solvers);
  set(gca, 'XTickLabel', {solvers.name});
  filename = sprintf('../figures/noise/%d-%g.eps', i, noise_stddev);
  print(filename, '-depsc2');
  unix(['epstopdf ', filename]);
  unix(['rm ', filename]);
end

residuals = permute(residuals, [1, 3, 2]);
