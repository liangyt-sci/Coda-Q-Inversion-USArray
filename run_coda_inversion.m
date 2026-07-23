function run_coda_inversion(cfg, mode)
%RUN_CODA_INVERSION Shared driver for the inversion scripts.

if ~exist(cfg.output_dir, 'dir')
    mkdir(cfg.output_dir);
end

reordered_file = fullfile(cfg.output_dir, 'reordered_data.txt');
event_file = fullfile(cfg.output_dir, 'event_index.txt');
station_file = fullfile(cfg.output_dir, 'station_index.txt');

data0 = reorder_coda_for_inversion(cfg.combined_data_file, reordered_file, ...
    cfg.n_amp, event_file, station_file);
data0 = data0(:, 1:(12 + cfg.n_amp));
data0 = data0(isfinite(data0(:,1)) & isfinite(data0(:,2)) & any(isfinite(data0(:,13:end)),2), :);
if isempty(data0)
    error('No valid rows remain after reordering.');
end

inv_cfg = cfg;
inv_cfg.mode = mode;

[data, model, diagnostics] = run_inversion_once(data0, inv_cfg);
write_outputs(cfg.output_dir, data, model, diagnostics, inv_cfg);
end

function [data, model, diagnostics] = run_inversion_once(data0, cfg)
fprintf('\n=== Inversion ===\n');
data0 = data0(:, 1:(12 + cfg.n_amp));
[data_filt, ~, ~] = filter_event_station_pairs(data0, ...
    cfg.min_stations_per_event, cfg.min_events_per_station);
if isempty(data_filt)
    error('No rows remain after event-station filtering.');
end

[data, model, diagnostics] = coda_q_inversion(data_filt, cfg);
diagnostics.iterations = 1;
diagnostics.negative_metric_rows = sum(data(:, end) < 0);
end

function write_outputs(output_dir, data, model, diagnostics, cfg)
base_col = 12 + cfg.n_amp;
eve = data(:,1);
sta = data(:,2);
tbeg = data(:,12);

save(fullfile(output_dir, 'final_data.txt'), 'data', '-ascii');
used_info = diagnostics.used_info; %#ok<NASGU>
save(fullfile(output_dir, 'used_records.txt'), 'used_info', '-ascii');

[record_result, station_result] = summarize_results(data, model, diagnostics, cfg);
save(fullfile(output_dir, 'record_results.txt'), 'record_result', '-ascii');
save(fullfile(output_dir, 'station_results.txt'), 'station_result', '-ascii');

fit_pair = build_pair_prediction(eve, sta, tbeg, model, cfg);
save(fullfile(output_dir, 'two_point_fit_vs_model.txt'), 'fit_pair', '-ascii');

fprintf('\nDone.\n');
fprintf('Final data: %s\n', fullfile(output_dir, 'final_data.txt'));
fprintf('Station results: %s\n', fullfile(output_dir, 'station_results.txt'));
fprintf('Records used: %d, events: %d, stations: %d\n', ...
    size(data,1), diagnostics.n_events, diagnostics.n_stations);
fprintf('Rows with err_sum < 0 retained as diagnostics: %d\n', diagnostics.negative_metric_rows);
fprintf('Energy columns: %d\n', base_col - 12);
end

function [record_result, station_result] = summarize_results(data, model, diagnostics, cfg)
eve = data(:,1);
sta = data(:,2);
stla = data(:,10);
stlo = data(:,11);
nrow = size(data,1);

param_std = estimate_parameter_std(model.A, norm(data(:,end-2)) / max(1, 2*nrow), cfg.uncertainty_rank);
evenum = diagnostics.n_events;
stanum = diagnostics.n_stations;
qc2_std = param_std(evenum+stanum+evenum+1 : end);

record_result = nan(nrow, 8);
record_result(:,1) = sta;
record_result(:,2) = stlo;
record_result(:,3) = stla;
record_result(:,4) = model.qc2_res(sta);
record_result(:,5) = model.sta_res(sta);
record_result(:,6) = qc2_std(sta);
record_result(:,7) = -model.x3_sum(:,1);
record_result(:,8) = model.qc1_res(eve);

unique_sta = unique(sta);
station_result = nan(numel(unique_sta), 8);
for i = 1:numel(unique_sta)
    idx = sta == unique_sta(i);
    station_result(i,1) = unique_sta(i);
    for j = 2:8
        station_result(i,j) = mean(record_result(idx,j), 'omitnan');
    end
end
end

function fit_pair = build_pair_prediction(eve, sta, tbeg, model, cfg)
C = log10(exp(1));
nrow = numel(eve);
fit_pair = nan(nrow*2, 4);
for i = 1:nrow
    if strcmpi(cfg.mode, 'fixed_window')
        tt = [tbeg(i); tbeg(i) + cfg.window_seconds];
    else
        tt = [tbeg(i); cfg.tmax_lapse];
    end
    rows = (i-1)*2 + (1:2);
    fit_pair(rows,1) = i;
    fit_pair(rows,2) = tt;
    fit_pair(rows,3) = model.b3(rows);
    fit_pair(rows,4) = model.eve_res(eve(i)) + model.sta_res(sta(i)) ...
        - C * tt .* (model.qc1_res(eve(i)) + model.qc2_res(sta(i)));
end
end

function param_std = estimate_parameter_std(A, variance_scale, rank0)
param_std = nan(size(A,2), 1);
k = min([rank0, min(size(A))-1]);
if k < 1
    return;
end
try
    [~, S, V] = svds(A, k);
    singvals = diag(S);
    valid = singvals > 0;
    Cm = variance_scale * V(:,valid) * diag(1 ./ singvals(valid).^2) * V(:,valid)';
    param_std = sqrt(diag(Cm));
catch ME
    warning('Uncertainty estimate skipped: %s', ME.message);
end
end
