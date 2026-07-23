function [data_out, event_names, station_names] = reorder_coda_for_inversion(infile, outfile, n_amp, event_file, station_file)
%REORDER_CODA_FOR_INVERSION Convert coda text output to numeric input.


if nargin < 3 || isempty(n_amp)
    error('n_amp is required.');
end

fid = fopen(infile, 'r');
if fid == -1
    error('Cannot open input file: %s', infile);
end
cleanup_obj = onCleanup(@() fclose(fid));

fmt = ['%s', repmat(' %f', 1, 6), ' %s', repmat(' %f', 1, 3), repmat(' %f', 1, n_amp)];
C = textscan(fid, fmt, 'CollectOutput', false);
clear cleanup_obj;

event_names = C{1};
evlo = C{2};
evla = C{3};
evdp = C{4};
dist = C{5};
az = C{6};
mag = C{7};
station_names_in = C{8};
stla = C{9};
stlo = C{10};
tbeg = C{11};

amp = nan(numel(event_names), n_amp);
for i = 1:n_amp
    amp(:,i) = C{11+i};
end

[event_names_unique, ~, event_idx] = unique(event_names, 'stable');
[station_names, ~, station_idx] = unique(station_names_in, 'stable');
event_numeric = str2double(event_names);
bad_event_id = ~isfinite(event_numeric);
event_numeric(bad_event_id) = event_idx(bad_event_id);

data_out = [event_idx, station_idx, event_numeric, evlo, evla, evdp, ...
    dist, az, mag, stla, stlo, tbeg, amp];

if nargin >= 2 && ~isempty(outfile)
    ensure_parent_dir(outfile);
    save(outfile, 'data_out', '-ascii');
end

if nargin >= 4 && ~isempty(event_file)
    write_text_lines(event_file, event_names_unique);
end
if nargin >= 5 && ~isempty(station_file)
    write_text_lines(station_file, station_names);
end
end

function write_text_lines(filename, lines)
ensure_parent_dir(filename);
fid = fopen(filename, 'w');
if fid == -1
    error('Cannot open file for writing: %s', filename);
end
cleanup_obj = onCleanup(@() fclose(fid));
for i = 1:numel(lines)
    fprintf(fid, '%s\n', lines{i});
end
clear cleanup_obj;
end

function ensure_parent_dir(filename)
parent = fileparts(filename);
if ~isempty(parent) && ~exist(parent, 'dir')
    mkdir(parent);
end
end
