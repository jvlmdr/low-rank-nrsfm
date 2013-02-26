num_sequences = size(residuals, 1);
num_solvers = numel(solvers);
num_ranks = numel(ranks);

residuals = permute(residuals, [1, 3, 2]);
residuals = reshape(residuals, [num_sequences, num_solvers * num_ranks]);

figure;
boxplot(residuals);
hold on;
plotSpread(residuals);
hold off;
grid on;

axis([0, num_solvers * num_ranks + 1, 0, max(residuals(:)) * 1.05]);
title('3D reconstruction error');
set(gca, 'XTick', 1:(num_solvers * num_ranks));
names = repmat({solvers(:).name}, [num_ranks, 1]);
set(gca, 'XTickLabel', names(:));
%filename = sprintf('../figures/noise/%d-%g.eps', i, noise_stddev);
%print(filename, '-depsc2');
%unix(['epstopdf ', filename]);
%unix(['rm ', filename]);

residuals = reshape(residuals, [num_sequences, num_solvers, num_ranks]);
residuals = permute(residuals, [1, 3, 2]);
