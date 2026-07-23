function [new_data, model, diagnostics] = coda_q_inversion(data, cfg)
%CODA_Q_INVERSION Two-point event/station coda-Q inversion.

cfg = fill_defaults(cfg);
C = log10(exp(1));
base_col = 12 + cfg.n_amp;

if size(data,2) < base_col
    error('Input data has %d columns, but at least %d are required.', size(data,2), base_col);
end

data_core = data(:,1:base_col);
eve = data_core(:,1);
sta = data_core(:,2);
tbeg = data_core(:,12);
amp = data_core(:,13:base_col);
nrow0 = size(data_core,1);

[keep, n_valid_pts, last_valid_time] = select_valid_records(tbeg, amp, cfg);
if ~any(keep)
    error('No records satisfy the inversion window criteria.');
end

data_keep = data_core(keep,:);
eve = eve(keep);
sta = sta(keep);
tbeg = tbeg(keep);
amp = amp(keep,:);
n_valid_pts = n_valid_pts(keep);
last_valid_time = last_valid_time(keep);
original_indices = find(keep);

nrow = size(data_keep,1);
ncol = 2;
b3 = zeros(nrow*ncol, 1);
x3_sum = nan(nrow, 2);
rmse_sum1 = nan(nrow, 1);
r_sum1 = nan(nrow, 1);
used_info = nan(nrow, 7);

for i = 1:nrow
    [t_vec, amp_vec] = record_vectors(tbeg(i), amp(i,:), cfg);
    b_local = amp_vec + cfg.alpha * log10(t_vec);
    A_local = [C * t_vec, ones(numel(t_vec), 1)];
    x3 = A_local \ b_local;
    x3_sum(i,:) = x3(:)';

    tt_pair = record_pair_times(tbeg(i), cfg);
    rows = (i-1)*ncol + (1:ncol);
    b3(rows) = x3(2) + x3(1) * C * tt_pair;

    residual1 = b_local - A_local*x3;
    rmse_sum1(i) = sqrt(sum(residual1.^2));
    denom = sum((b_local - mean(b_local)).^2);
    if denom > 0
        r_sum1(i) = 1 - sum(residual1.^2) / denom;
    end

    used_info(i,:) = [original_indices(i), eve(i), sta(i), tbeg(i), ...
        tt_pair(end), n_valid_pts(i), last_valid_time(i)];
end

evenum = max(eve);
stanum = max(sta);
A = build_global_matrix(eve, sta, tbeg, cfg, evenum, stanum);
[x, flag, relres] = lsqr(A, b3, cfg.lsqr_tol, cfg.lsqr_maxit);

model.eve_res = x(1:evenum);
model.sta_res = x(evenum+1:evenum+stanum);
model.qc1_res = x(evenum+stanum+1:evenum+stanum+evenum);
model.qc2_res = x(evenum+stanum+evenum+1:end);
model.A = A;
model.x = x;
model.x3_sum = x3_sum;
model.b3 = b3;

rmse_sum2 = nan(nrow, 1);
for i = 1:nrow
    [t_vec, amp_vec] = record_vectors(tbeg(i), amp(i,:), cfg);
    b_obs = amp_vec + cfg.alpha * log10(t_vec);
    b_pred = model.eve_res(eve(i)) + model.sta_res(sta(i)) ...
        - C * t_vec .* (model.qc1_res(eve(i)) + model.qc2_res(sta(i)));
    rmse_sum2(i) = norm(b_obs - b_pred);
end

err_sum = rmse_sum1 * cfg.err_weight - rmse_sum2;
new_data = [data_keep, rmse_sum1, rmse_sum2, err_sum];

diagnostics.cfg = cfg;
diagnostics.used_info = used_info;
diagnostics.r_sum1 = r_sum1;
diagnostics.lsqr_flag = flag;
diagnostics.lsqr_relres = relres;
diagnostics.n_original = nrow0;
diagnostics.n_used = nrow;
diagnostics.n_events = evenum;
diagnostics.n_stations = stanum;
end

function cfg = fill_defaults(cfg)
if ~isfield(cfg, 'mode') || isempty(cfg.mode), cfg.mode = 'fixed_window'; end
if ~isfield(cfg, 'n_amp') || isempty(cfg.n_amp), error('cfg.n_amp is required.'); end
if ~isfield(cfg, 'dt_energy') || isempty(cfg.dt_energy), cfg.dt_energy = 0.5; end
if ~isfield(cfg, 'window_seconds') || isempty(cfg.window_seconds), cfg.window_seconds = 40; end
if ~isfield(cfg, 'tmax_lapse') || isempty(cfg.tmax_lapse), cfg.tmax_lapse = 100; end
if ~isfield(cfg, 'min_coda_seconds') || isempty(cfg.min_coda_seconds), cfg.min_coda_seconds = 30; end
if ~isfield(cfg, 'alpha') || isempty(cfg.alpha), cfg.alpha = 1.5; end
if ~isfield(cfg, 'err_weight') || isempty(cfg.err_weight), cfg.err_weight = 60; end
if ~isfield(cfg, 'lsqr_tol') || isempty(cfg.lsqr_tol), cfg.lsqr_tol = 1e-8; end
if ~isfield(cfg, 'lsqr_maxit') || isempty(cfg.lsqr_maxit), cfg.lsqr_maxit = 2000; end
end

function [keep, n_valid_pts, last_valid_time] = select_valid_records(tbeg, amp, cfg)
nrow = size(amp,1);
keep = false(nrow,1);
n_valid_pts = zeros(nrow,1);
last_valid_time = nan(nrow,1);
min_pts = ceil(cfg.min_coda_seconds / cfg.dt_energy);

for i = 1:nrow
    [t_vec, amp_vec] = record_vectors(tbeg(i), amp(i,:), cfg);
    if numel(t_vec) < min_pts
        continue;
    end
    if max(t_vec) - min(t_vec) < cfg.min_coda_seconds
        continue;
    end
    if any(~isfinite(amp_vec))
        continue;
    end
    keep(i) = true;
    n_valid_pts(i) = numel(t_vec);
    last_valid_time(i) = max(t_vec);
end
end

function [t_vec, amp_vec] = record_vectors(tbeg, amp_row, cfg)
t_full = tbeg + (0:cfg.n_amp-1) * cfg.dt_energy;
switch lower(cfg.mode)
    case 'fixed_window'
        t_end = tbeg + cfg.window_seconds;
    case {'fixed_tmax', 'fix_max', 'fixed_max'}
        t_end = cfg.tmax_lapse;
    otherwise
        error('Unknown inversion mode: %s', cfg.mode);
end
valid = isfinite(t_full) & t_full > 0 & t_full <= t_end & isfinite(amp_row);
t_vec = t_full(valid)';
amp_vec = amp_row(valid)';
end

function tt_pair = record_pair_times(tbeg, cfg)
switch lower(cfg.mode)
    case 'fixed_window'
        tt_pair = [tbeg; tbeg + cfg.window_seconds];
    case {'fixed_tmax', 'fix_max', 'fixed_max'}
        tt_pair = [tbeg; cfg.tmax_lapse];
end
end

function A = build_global_matrix(eve, sta, tbeg, cfg, evenum, stanum)
C = log10(exp(1));
nrow = numel(eve);
ncol = 2;
b_nlen = nrow * ncol;
a_ncol = evenum + stanum + evenum + stanum;

I = zeros(nrow*ncol*4, 1);
J = zeros(nrow*ncol*4, 1);
V = zeros(nrow*ncol*4, 1);
nz = 0;

for i = 1:nrow
    rows = (i-1)*ncol + (1:ncol);
    tt_pair = record_pair_times(tbeg(i), cfg);

    for k = 1:ncol
        row_id = rows(k);
        t_now = tt_pair(k);

        nz = nz + 1; I(nz) = row_id; J(nz) = eve(i); V(nz) = 1;
        nz = nz + 1; I(nz) = row_id; J(nz) = evenum + sta(i); V(nz) = 1;
        nz = nz + 1; I(nz) = row_id; J(nz) = evenum + stanum + eve(i); V(nz) = -C * t_now;
        nz = nz + 1; I(nz) = row_id; J(nz) = evenum + stanum + evenum + sta(i); V(nz) = -C * t_now;
    end
end

A = sparse(I(1:nz), J(1:nz), V(1:nz), b_nlen, a_ncol);
end
