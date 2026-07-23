function [new_data, kept_events, kept_stations] = filter_event_station_pairs(data, min_stations_per_event, min_events_per_station)
%FILTER_EVENT_STATION_PAIRS Keep a connected event-station subset.
%
% data(:,1) must be event indices, and data(:,2) station indices. The
% function iteratively removes events and stations that do not meet the
% requested minimum coverage, then renumbers the retained event/station
% indices to consecutive values.

if nargin < 2 || isempty(min_stations_per_event)
    min_stations_per_event = 5;
end
if nargin < 3 || isempty(min_events_per_station)
    min_events_per_station = 5;
end

event_ids = data(:,1);
station_ids = data(:,2);
keep_data = true(size(data,1),1);
changed = true;

while changed
    changed = false;

    [unique_events, ~, event_idx] = unique(event_ids(keep_data));
    event_counts = accumarray(event_idx, 1);
    bad_events = unique_events(event_counts < min_stations_per_event);

    [unique_stations, ~, station_idx] = unique(station_ids(keep_data));
    station_counts = accumarray(station_idx, 1);
    bad_stations = unique_stations(station_counts < min_events_per_station);

    remove_now = keep_data & (ismember(event_ids, bad_events) | ismember(station_ids, bad_stations));
    if any(remove_now)
        keep_data(remove_now) = false;
        changed = true;
    end
end

filtered_data = data(keep_data,:);
if isempty(filtered_data)
    new_data = filtered_data;
    kept_events = [];
    kept_stations = [];
    return;
end

kept_events = unique(filtered_data(:,1));
kept_stations = unique(filtered_data(:,2));

[~, new_eve] = ismember(filtered_data(:,1), kept_events);
[~, new_sta] = ismember(filtered_data(:,2), kept_stations);
new_data = [new_eve, new_sta, filtered_data(:,3:end)];

fprintf('Event-station filter: %d -> %d records, %d events, %d stations.\n', ...
    size(data,1), size(new_data,1), numel(kept_events), numel(kept_stations));
end
