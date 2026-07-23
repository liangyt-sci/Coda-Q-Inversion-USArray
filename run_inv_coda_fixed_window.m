%% Fixed-window coda-Q inversion
% Edit the configuration block, then run this script.

clear; clc;

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);

cfg = struct();
cfg.combined_data_file = '/path/to/combined_data.txt';
cfg.output_dir = '/path/to/inversion_fixed_window_output';

cfg.window_seconds = 40;
cfg.dt_energy = 0.5;
cfg.n_amp = round(cfg.window_seconds / cfg.dt_energy);
cfg.min_coda_seconds = 30;

cfg.min_stations_per_event = 5;
cfg.min_events_per_station = 5;

cfg.alpha = 1.5;
cfg.err_weight = 60;
cfg.uncertainty_rank = 1000;

run_coda_inversion(cfg, 'fixed_window');
