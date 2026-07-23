%% coda waveform processing
% Measure coda-window mean-square energy from SAC waveforms.


clc;

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);

if exist('CODA_CFG', 'var')
    cfg = CODA_CFG;
    clear CODA_CFG;
else
    cfg = struct();
    cfg.waveform_root = '/path/to/waveforms';     % event folders live here
    cfg.event_list    = '/path/to/eve_list.txt';  % one event id/folder per line
    cfg.station_list  = '/path/to/sta_list.txt';  % one station code per line
    cfg.output_root   = '/path/to/coda_output';
    cfg.figure_root   = '/path/to/coda_figures';

    cfg.sampling_rate = 40.0;
    cfg.taper_seconds = 2.0;
    cfg.smooth_ms_nt = 160;   % samples in each moving mean-square window
    cfg.smooth_nt_int = 20;   % samples between adjacent energy points
    cfg.active_frequency_labels = {'1.5Hz','3Hz','6Hz','12Hz'};

    cfg.regions = define_regions();
    cfg.active_region_names = {'all'};

    cfg.coda_window_seconds = 40;
    cfg.fixed_tmax_seconds = [];  % [] means fixed length; set e.g. 100 for fixed Tmax/fix max
    cfg.min_coda_seconds = 30;    % used only when fixed_tmax_seconds is not empty
    cfg.noise_window_seconds = 5;
    cfg.origin_time_seconds = 30;
    cfg.min_snr = 3;
    cfg.min_fit_corr = 0.9;
    cfg.max_event_depth_km = 40;
    cfg.pad_fixed_tmax = true;    % fixed Tmax output has a constant number of energy columns
    cfg.save_figures = false;
    cfg.verbose = false;
end

run_coda_processing(cfg);

function run_coda_processing(cfg)
old_visibility = get(groot, 'DefaultFigureVisible');
cleanup_obj = onCleanup(@() set(groot, 'DefaultFigureVisible', old_visibility));
if ~cfg.save_figures
    set(groot, 'DefaultFigureVisible', 'off');
end

frequencies = select_by_name(define_frequency_configs(cfg), ...
    cfg.active_frequency_labels, 'frequency');
regions = select_by_name(cfg.regions, cfg.active_region_names, 'region');
stations = read_text_lines(cfg.station_list);
events = read_text_lines(cfg.event_list);

for ifreq = 1:numel(frequencies)
    freq = frequencies(ifreq);
    [b, a] = butter(4, freq.band_hz / (cfg.sampling_rate / 2));

    for iregion = 1:numel(regions)
        region = regions(iregion);
        if isempty(cfg.fixed_tmax_seconds)
            mode_name = 'fixed_window';
        else
            mode_name = sprintf('fixed_tmax_%gs', cfg.fixed_tmax_seconds);
        end
    outdir = fullfile(cfg.output_root, mode_name, ['coda_', freq.label], region.name);
        figdir = fullfile(cfg.figure_root, freq.label, region.name);
        ensure_dir(outdir);
        if cfg.save_figures
            ensure_dir(figdir);
        end

        fprintf('Processing %s / %s\n', freq.label, region.name);
        for ista = 1:numel(stations)
            station = stations{ista};
            outfile = fullfile(outdir, [station, '.txt']);
            fid = fopen(outfile, 'w');
            if fid == -1
                error('Cannot open output file: %s', outfile);
            end
            close_file = onCleanup(@() fclose(fid));

            for ievent = 1:numel(events)
                event_id = events{ievent};
                event_dir = fullfile(cfg.waveform_root, event_id);
                sac_files = dir(fullfile(event_dir, [station, '*BH*.sac']));
                if isempty(sac_files)
                    continue;
                end

                result = process_one_record(event_dir, sac_files, station, event_id, ...
                    b, a, freq, region, cfg);
                if ~result.keep
                    if cfg.verbose && ~isempty(result.reason)
                        fprintf('Skip %s %s: %s\n', event_id, station, result.reason);
                    end
                    continue;
                end

                fprintf(fid, '%s %8.4f %8.4f %8.4f %8.3f %8.3f %5.2f %s %8.4f %8.4f %8.4f ', ...
                    event_id, result.evlo, result.evla, result.evdp, result.dist_km, ...
                    result.az, result.mag, station, result.stla, result.stlo, result.tbeg);
                fprintf(fid, '%10.5f ', result.log_energy_out);
                fprintf(fid, '\n');

                if cfg.save_figures
                    plot_coda_trace(result.trace, figdir, cfg.save_figures);
                end
            end
            clear close_file;
        end
    end
end
end

function result = process_one_record(event_dir, sac_files, station, event_id, b, a, freq, region, cfg)
result = struct('keep', false, 'reason', '', 'energy', [], 'trace', []);
s = readsac(fullfile(event_dir, sac_files(1).name));

if s.NPTS < 500 || isnan(s.NPTS)
    result.reason = 'too few samples';
    return;
end
if ~in_region(s, region)
    result.reason = 'outside region';
    return;
end
if s.T1 < 0 || s.T2 < 0 || s.T1 < cfg.origin_time_seconds || s.E < 200
    result.reason = 'invalid picks or short record';
    return;
end

delta = s.DELTA;
raw = detrend(s.DATA1(:));
raw = hann_taper(raw, delta, s.NPTS, cfg.taper_seconds);
filtered = filter(b, a, raw);

nnoise = floor(cfg.noise_window_seconds / delta);
noise_beg = floor(s.T1 / delta) - nnoise;
noise_end = floor(s.T1 / delta) - 1;
if noise_beg < 1 || noise_end > numel(filtered) || noise_end <= noise_beg
    result.reason = 'noise window out of range';
    return;
end
noise_level = mean(filtered(noise_beg:noise_end).^2);

s.O = cfg.origin_time_seconds;
tbeg_abs = s.T2 + s.T2 - s.O;
tbeg_lapse = tbeg_abs - s.O;
nbeg = floor((tbeg_abs - s.B) / delta) + 1;
if nbeg < 1 || ~isfinite(nbeg)
    result.reason = 'coda start out of range';
    return;
end

if isempty(cfg.fixed_tmax_seconds)
    s_nt = floor(cfg.coda_window_seconds / (freq.nt_int * delta) + 0.5);
    if ~smoothing_window_in_range(nbeg, s_nt, freq.ms_nt, freq.nt_int, s.NPTS)
        result.reason = 'coda smoothing window out of range';
        return;
    end
    energy_raw = mean_square_coda(filtered, nbeg, s_nt, freq.ms_nt, freq.nt_int);
    time_energy = tbeg_abs + (0:s_nt-1)' * freq.nt_int * delta;
else
    dt_energy = freq.nt_int * delta;
    coda_len = cfg.fixed_tmax_seconds - tbeg_lapse;
    if ~isfinite(coda_len) || coda_len < cfg.min_coda_seconds
        result.reason = 'fixed Tmax coda window too short';
        return;
    end
    s_nt = floor(coda_len / dt_energy);
    s_nt_max = ceil(cfg.fixed_tmax_seconds / dt_energy);
    if s_nt < ceil(cfg.min_coda_seconds / dt_energy)
        result.reason = 'too few fixed Tmax energy points';
        return;
    end
    if ~smoothing_window_in_range(nbeg, s_nt, freq.ms_nt, freq.nt_int, s.NPTS)
        result.reason = 'fixed Tmax smoothing window out of range';
        return;
    end
    [energy_raw, time_energy, ok] = mean_square_coda_fixed_tmax(filtered, nbeg, tbeg_lapse, ...
        s_nt, freq.ms_nt, freq.nt_int, delta, cfg.fixed_tmax_seconds);
    if ~ok
        result.reason = 'fixed Tmax energy extraction failed';
        return;
    end
    energy_tmp = energy_raw - noise_level;
    valid_energy = isfinite(energy_tmp) & energy_tmp > 0;
    if sum(valid_energy) / numel(energy_tmp) < 0.7
        result.reason = 'too many non-positive fixed Tmax energy values';
        return;
    end
    log_energy_out = nan(s_nt_max, 1);
    log_energy_out(1:numel(energy_tmp)) = log10(max(energy_tmp, eps));
end

snr = mean(energy_raw) / noise_level;
energy = energy_raw - noise_level;
if snr < cfg.min_snr || any(energy <= 0) || s.EVDP >= cfg.max_event_depth_km
    result.reason = 'failed snr/energy/depth criteria';
    return;
end
if isempty(cfg.fixed_tmax_seconds)
    log_energy_out = log10(energy);
end

x = (1:numel(energy))';
fit_coef = polyfit(x, log10(energy), 1);
fit_log = fit_coef(1) * x + fit_coef(2);
fit_corr = corr(fit_log, log10(energy));
if fit_coef(1) > 0 && abs(fit_corr) < cfg.min_fit_corr
    result.reason = 'positive weak coda decay';
    return;
end

[dist_km, az] = great_circle_distance_km(s.EVLA, s.EVLO, s.STLA, s.STLO);
if ~isfinite(dist_km)
    result.reason = 'invalid distance';
    return;
end

result.keep = true;
result.energy = energy(:);
result.log_energy_out = log_energy_out(:);
result.evlo = s.EVLO;
result.evla = s.EVLA;
result.evdp = s.EVDP;
result.dist_km = dist_km;
result.az = az;
result.mag = 0;
result.stla = s.STLA;
result.stlo = s.STLO;
if isempty(cfg.fixed_tmax_seconds)
    result.tbeg = tbeg_abs;
else
    result.tbeg = tbeg_lapse;
end

components = read_components(event_dir, sac_files, b, a, cfg.taper_seconds);
result.trace = build_trace_struct(s, event_id, station, raw, filtered, components, ...
    energy_raw, energy, time_energy, fit_log, fit_corr, snr, tbeg_abs, cfg, freq, ...
    noise_beg, noise_end, result.dist_km);
end

function arr_ms = mean_square_coda(arrin, nbeg, s_nt, ms_nt, nt_int)
arr_sq = arrin(:).^2;
arr_ms = zeros(s_nt, 1);
half_win = floor(ms_nt / 2);
for i = 1:s_nt
    i1 = nbeg + (i-1) * nt_int - half_win;
    i2 = i1 + ms_nt - 1;
    arr_ms(i) = mean(arr_sq(i1:i2));
end
end

function [arr_ms, t_ms, ok] = mean_square_coda_fixed_tmax(arrin, nbeg, tbeg, s_nt, ms_nt, nt_int, delta, tmax)
ok = true;
arr_sq = arrin(:).^2;
arr_ms_tmp = nan(s_nt, 1);
t_ms_tmp = nan(s_nt, 1);
half_win = floor(ms_nt / 2);
end_offset = ms_nt - 1 - half_win;
n_valid = 0;
for i = 1:s_nt
    ncenter = nbeg + (i-1) * nt_int;
    t_center = tbeg + (i-1) * nt_int * delta;
    if t_center + end_offset * delta > tmax
        break;
    end
    i1 = ncenter - half_win;
    i2 = i1 + ms_nt - 1;
    if i1 < 1 || i2 > numel(arr_sq)
        ok = false;
        arr_ms = [];
        t_ms = [];
        return;
    end
    n_valid = n_valid + 1;
    arr_ms_tmp(n_valid) = mean(arr_sq(i1:i2));
    t_ms_tmp(n_valid) = t_center;
end
ok = n_valid > 0;
arr_ms = arr_ms_tmp(1:n_valid);
t_ms = t_ms_tmp(1:n_valid);
end

function ok = smoothing_window_in_range(nbeg, s_nt, ms_nt, nt_int, npts)
half_win = floor(ms_nt / 2);
first_needed = nbeg - half_win;
last_needed = nbeg + (s_nt - 1) * nt_int + (ms_nt - 1 - half_win);
ok = first_needed >= 1 && last_needed <= npts;
end

function trace = build_trace_struct(s, event_id, station, raw, filtered, components, ...
    energy_raw, energy, time_energy, fit_log, fit_corr, snr, tbeg, cfg, freq, ...
    noise_beg, noise_end, dist_km)
trace = struct();
trace.event_id = strtrim(event_id);
trace.station = strtrim(station);
trace.freq_label = freq.label;
trace.band_hz = freq.band_hz;
trace.t = s.B + (0:s.NPTS-1)' * s.DELTA;
trace.raw = raw(:);
trace.filtered = filtered(:);
trace.components = components;
trace.energy_log = log10(energy(:));
trace.energy_raw_log = log10(max(energy_raw(:), eps));
trace.energy_time = time_energy(:);
trace.fit_log = fit_log(:);
trace.fit_corr = fit_corr;
trace.snr = snr;
trace.origin = s.O;
trace.t1 = s.T1;
trace.t2 = s.T2;
trace.tbeg = tbeg;
if isempty(cfg.fixed_tmax_seconds)
    trace.tend = time_energy(end);
else
    trace.tend = s.O + time_energy(end);
end
trace.noise_beg = s.B + (noise_beg - 1) * s.DELTA;
trace.noise_end = s.B + (noise_end - 1) * s.DELTA;
trace.dist_km = dist_km;
trace.window_seconds = cfg.coda_window_seconds;
end

function components = read_components(event_dir, sac_files, b, a, taper_seconds)
order = 'ZNE';
components = repmat(struct('name', '', 't', [], 'data', []), 3, 1);
for i = 1:3
    components(i).name = order(i);
end
for i = 1:numel(sac_files)
    token = regexp(sac_files(i).name, 'BH([ZNE])', 'tokens', 'once');
    if isempty(token)
        continue;
    end
    idx = find(order == token{1}, 1);
    sc = readsac(fullfile(event_dir, sac_files(i).name));
    data = detrend(sc.DATA1(:));
    data = hann_taper(data, sc.DELTA, sc.NPTS, taper_seconds);
    components(idx).t = sc.B + (0:sc.NPTS-1)' * sc.DELTA;
    components(idx).data = filter(b, a, data);
end
first = find(~cellfun(@isempty, {components.data}), 1);
if isempty(first)
    return;
end
for i = 1:3
    if isempty(components(i).data)
        components(i).t = components(first).t;
        components(i).data = zeros(size(components(first).data));
    end
end
end

function [dist_km, azimuth_deg] = great_circle_distance_km(lat1, lon1, lat2, lon2)
if any(~isfinite([lat1, lon1, lat2, lon2]))
    dist_km = NaN;
    azimuth_deg = NaN;
    return;
end
earth_radius_km = 6371.0088;
phi1 = deg2rad(lat1);
phi2 = deg2rad(lat2);
dphi = deg2rad(lat2 - lat1);
dlambda = deg2rad(lon2 - lon1);
a = sin(dphi/2).^2 + cos(phi1).*cos(phi2).*sin(dlambda/2).^2;
c = 2 * atan2(sqrt(a), sqrt(max(0, 1-a)));
dist_km = earth_radius_km * c;

y = sin(dlambda) * cos(phi2);
x = cos(phi1) * sin(phi2) - sin(phi1) * cos(phi2) * cos(dlambda);
azimuth_deg = mod(rad2deg(atan2(y, x)) + 360, 360);
end

function freqs = define_frequency_configs(cfg)
freqs = struct('name', {}, 'label', {}, 'center_hz', {}, 'band_hz', {}, 'ms_nt', {}, 'nt_int', {});
freqs(end+1) = make_frequency('1.5Hz', 1.5, [1 2], cfg);
freqs(end+1) = make_frequency('3Hz', 3.0, [2 4], cfg);
freqs(end+1) = make_frequency('6Hz', 6.0, [4 8], cfg);
freqs(end+1) = make_frequency('12Hz', 12.0, [8 16], cfg);
end

function freq = make_frequency(label, center_hz, band_hz, cfg)
freq.name = label;
freq.label = label;
freq.center_hz = center_hz;
freq.band_hz = band_hz;
freq.ms_nt = cfg.smooth_ms_nt;
freq.nt_int = cfg.smooth_nt_int;
end

function regions = define_regions()
regions = struct('name', {}, 'event_lon', {}, 'event_lat', {}, 'station_lon', {}, 'station_lat', {});
regions(end+1) = make_region('all', [], [], [], []);

end

function region = make_region(name, event_lon, event_lat, station_lon, station_lat)
region.name = name;
region.event_lon = event_lon;
region.event_lat = event_lat;
region.station_lon = station_lon;
region.station_lat = station_lat;
end

function tf = in_region(s, region)
tf = in_box(s.EVLO, s.EVLA, region.event_lon, region.event_lat) && ...
     in_box(s.STLO, s.STLA, region.station_lon, region.station_lat);
end

function tf = in_box(lon, lat, lon_lim, lat_lim)
if isempty(lon_lim) || isempty(lat_lim)
    tf = true;
else
    tf = isfinite(lon) && isfinite(lat) && lon >= lon_lim(1) && lon <= lon_lim(2) && ...
        lat >= lat_lim(1) && lat <= lat_lim(2);
end
end

function selected = select_by_name(items, names, label)
if isempty(names)
    selected = items;
    return;
end
selected = items([]);
for i = 1:numel(names)
    idx = find(strcmp({items.name}, names{i}), 1);
    if isempty(idx)
        error('Unknown %s: %s', label, names{i});
    end
    selected(end+1) = items(idx);
end
end

function lines = read_text_lines(filename)
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file: %s', filename);
end
cleanup_obj = onCleanup(@() fclose(fid));
lines = {};
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if ischar(line) && ~isempty(line)
        lines{end+1, 1} = line; 
    end
end
clear cleanup_obj;
end

function ensure_dir(dirname)
if ~exist(dirname, 'dir')
    mkdir(dirname);
end
end
