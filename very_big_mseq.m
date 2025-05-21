clear; clc;
close all;

c = physconst("Lightspeed");

targets = struct( ...
    'elevation', [15], ...
    'azimuth', [25], ...
    'range', [1000], ...
    'speed', []);

max_time = max(targets.range) * 2 /c;

fc = 24e9;
lambda = c/fc;
d = lambda;

%% Tx-Rx array
Ntx = 8; Nrx = 8; % TODO: сейчас количество передатчиков ровно чтобы М 
                % последовательность делилась на числа передатчиков
tx_e = (0 : (Ntx - 1)) * d;
rx_e = ((0 : (Nrx - 1)) * d)';

%% Generate M-seq
n = 12;
ps_phases = get_mseq_n_times(n, 1);
pack_len = fix((2^n-1)/Ntx);
new_pss = zeros(Ntx, pack_len);
for i = 1 : Ntx - 1
    start_index = ((i-1) * pack_len) + 1;
    end_index = start_index + pack_len - 1;
    new_pss(i, :) = ps_phases(start_index : end_index);
end

abs_scan_limit = 1/lambda;
scanning_theta = -45 : 1 : 50;
scanning_phi = -34 : 1 : 40;

beams = mimo_model(new_pss, tx_e, rx_e, scanning_phi, scanning_theta, targets, lambda);

%% Find angles of maximums
tmp_plane_data = max(abs(beams), [], 3); % required to find azimuth and elevation of maximums
[est_azim, est_elev, azim_idx, elev_idx] = find_max_direction(tmp_plane_data, scanning_phi, scanning_theta)
clear tmp_plane_data;
%% find range maximum
[r_max, r_idx] = max(beams(elev_idx, azim_idx, :));

%% filter maximum
window = 2;

filtered_beams = beams;
filtered_beams(:, :, r_idx - window / 2 : r_idx + window / 2) = 0;

%% Results
all_max = max(abs(beams), [], [3, 2, 1])
all_rms = rms(beams, [3, 2, 1])
all_mean = mean(abs(beams), [3, 2, 1])

all_filtered_max = max(abs(filtered_beams), [], [3, 2, 1])
all_filtered_rms = rms(filtered_beams, [3, 2, 1])
all_filtered_mean = mean(abs(filtered_beams), [3, 2, 1])