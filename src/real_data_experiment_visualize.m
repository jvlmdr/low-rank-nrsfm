num_sequences = numel(scenes);
num_ranks = numel(ranks);
num_cameras_sans = numel(solvers.find_cameras);
num_cameras_with = numel(solvers.find_cameras_with_K);
num_structure_sans = numel(solvers.find_structure);
num_structure_with = numel(solvers.find_structure_with_K);
num_frames = size(errors.sans_sans, 1);

mean_errors = struct(...
    'sans_sans', shiftdim(mean(errors.sans_sans), 1), ...
    'sans_with', shiftdim(mean(errors.sans_with), 1), ...
    'with_sans', shiftdim(mean(errors.with_sans), 1), ...
    'with_with', shiftdim(mean(errors.with_with), 1));

mean_benchmarks = shiftdim(mean(benchmarks));

best_errors = struct(...
    'sans_sans', mean_errors.sans_sans, ...
    'sans_with', shiftdim(min(mean_errors.sans_with), 1), ...
    'with_sans', shiftdim(min(mean_errors.with_sans), 1), ...
    'with_with', shiftdim(min(mean_errors.with_with), 1));

all_errors = [...
  reshape(best_errors.sans_sans, [], num_sequences); ...
  reshape(best_errors.sans_with, [], num_sequences); ...
  reshape(best_errors.with_sans, [], num_sequences); ...
  reshape(best_errors.with_with, [], num_sequences); ...
]';

n = num_cameras_sans * num_structure_sans + ...
    num_cameras_sans * num_structure_with + ...
    num_cameras_with * num_structure_sans + ...
    num_cameras_with * num_structure_with;

all_errors = bsxfun(@rdivide, all_errors, mean_benchmarks);
all_errors = log10(all_errors);

figure;
boxplot(all_errors);
hold on;
plotSpread(all_errors);
hold off;
grid on;

ymin = min(all_errors(:));
ymax = max(all_errors(:));
ylim = [ymin - (ymax - ymin) * 0.1, ymax + (ymax - ymin) * 0.1];
axis([0, n + 1, ylim]);

%title('3D reconstruction error');
%set(gca, 'XTick', 1:(num_solvers * num_ranks));
%names = repmat({solvers(:).name}, [num_ranks, 1]);
%set(gca, 'XTickLabel', names(:));
%%filename = sprintf('../figures/noise/%d-%g.eps', i, noise_stddev);
%%print(filename, '-depsc2');
%%unix(['epstopdf ', filename]);
%%unix(['rm ', filename]);
