function plot_coda_trace(trace, outdir, keep_open)
%PLOT_CODA_TRACE Save a two-panel waveform and coda-fit quality-control plot.

if nargin < 3 || isempty(keep_open)
    keep_open = false;
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

fig = figure('Color', 'w', 'Position', [80 80 1650 560], ...
    'Visible', ternary(keep_open, 'on', 'off'));
left_ax = axes('Parent', fig, 'Position', [0.06 0.14 0.43 0.74]);
right_ax = axes('Parent', fig, 'Position', [0.56 0.14 0.38 0.74]);

plot_wave_panel(left_ax, trace);
plot_energy_panel(right_ax, trace);

annotation(fig, 'textbox', [0.08 0.93 0.84 0.045], ...
    'String', sprintf('%s | Event %s | Station %s | %.1f km', ...
    trace.freq_label, trace.event_id, trace.station, trace.dist_km), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontWeight', 'bold', 'FontSize', 13, 'EdgeColor', 'none');

safe_event = regexprep(trace.event_id, '[^\w.-]', '_');
safe_sta = regexprep(trace.station, '[^\w.-]', '_');
outfile = fullfile(outdir, sprintf('%s_%s_%s.png', safe_event, safe_sta, trace.freq_label));
print(fig, outfile, '-dpng', '-r300', '-noui');

if keep_open
    drawnow;
else
    close(fig, 'force');
end
end

function plot_wave_panel(ax, trace)
hold(ax, 'on');
time_zero = trace.origin;
plot_start = max(trace.origin - 10, trace.t(1));
plot_end = min(trace.tend + 25, trace.t(end));
yl = [0.35 4.65];

patch(ax, [trace.tbeg trace.tend trace.tend trace.tbeg] - time_zero, ...
    [yl(1) yl(1) yl(2) yl(2)], [0.98 0.90 0.70], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.45);
patch(ax, [trace.noise_beg trace.noise_end trace.noise_end trace.noise_beg] - time_zero, ...
    [3.00 3.00 3.42 3.42], [0.55 0.65 1.00], ...
    'EdgeColor', [0.20 0.30 0.95], 'FaceAlpha', 0.12);

energy_idx = trace.t >= plot_start & trace.t <= plot_end;
if any(energy_idx)
    raw_energy = trace.filtered(energy_idx).^2;
    raw_energy = raw_energy / max(max(raw_energy), eps);
    plot(ax, trace.t(energy_idx) - time_zero, raw_energy * 0.45 + 4.05, ...
        'Color', [0.20 0.20 0.20], 'LineWidth', 1.0);
end

offsets = [3.2 2.2 1.2];
labels = {'Z component', 'N component', 'E component'};
for i = 1:3
    tr = trace.components(i);
    idx = tr.t >= plot_start & tr.t <= plot_end;
    if any(idx)
        scale = max(abs(tr.data(idx)));
        plot(ax, tr.t(idx) - time_zero, tr.data(idx) / max(scale, eps) * 0.42 + offsets(i), ...
            'Color', [0.78 0.80 0.82], 'LineWidth', 0.8);
    end
    text(ax, plot_end - time_zero - 0.06 * (plot_end - plot_start), offsets(i) + 0.18, labels{i}, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'FontSize', 9);
end

plot(ax, [trace.origin trace.origin] - time_zero, yl, 'k-', 'LineWidth', 0.9);
plot(ax, [trace.t1 trace.t1] - time_zero, yl, '-', 'Color', [0.35 0.55 0.95], 'LineWidth', 1.1);
plot(ax, [trace.t2 trace.t2] - time_zero, yl, '-', 'Color', [0.88 0.48 0.14], 'LineWidth', 1.1);
plot(ax, [trace.tbeg trace.tbeg] - time_zero, yl, '--', 'Color', [0.55 0.55 0.55]);
plot(ax, [trace.tend trace.tend] - time_zero, yl, '--', 'Color', [0.55 0.55 0.55]);

xlim(ax, [plot_start plot_end] - time_zero);
ylim(ax, yl);
set(ax, 'YTick', []);
xlabel(ax, 'Time from origin (s)');
title(ax, 'Waveform and coda window');
grid(ax, 'on');
box(ax, 'on');
hold(ax, 'off');
end

function plot_energy_panel(ax, trace)
hold(ax, 'on');
plot(ax, trace.energy_time, trace.energy_raw_log, '-', ...
    'Color', [1.00 0.55 0.48], 'LineWidth', 0.9, 'DisplayName', 'raw smoothed energy');
plot(ax, trace.energy_time, trace.energy_log, 'o-', ...
    'Color', [0.16 0.38 0.95], 'MarkerFaceColor', 'w', 'MarkerSize', 4.0, ...
    'LineWidth', 1.4, 'DisplayName', 'noise-corrected energy');
plot(ax, trace.energy_time, trace.fit_log, '-', ...
    'Color', [0.28 0.78 0.24], 'LineWidth', 2.0, 'DisplayName', 'linear fit');

xlabel(ax, 'Lapse time (s)');
ylabel(ax, 'log_{10}(energy)');
title(ax, sprintf('r = %.2f, SNR = %.1f', trace.fit_corr, trace.snr));
legend(ax, 'Location', 'northwest', 'Box', 'off');
grid(ax, 'on');
box(ax, 'on');
hold(ax, 'off');
end

function out = ternary(cond, true_value, false_value)
if cond
    out = true_value;
else
    out = false_value;
end
end
