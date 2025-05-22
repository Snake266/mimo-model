function [beams, filtered_beams, all_max, all_filtered_max, all_filtered_rms] = mimo_system(targets, tx_signals, Ntx, Nrx, scanning_phi, scanning_theta)

c = physconst("Lightspeed");

fc = 24e9;
lambda = c/fc;
d = lambda;

%% Tx-Rx array
tx_e = (0 : (Ntx - 1)) * d;
rx_e = ((0 : (Nrx - 1)) * d)';

beams = mimo_model(tx_signals, tx_e, rx_e, scanning_phi, scanning_theta, targets, lambda);

%% Find angles of maximums
tmp_plane_data = max(abs(beams), [], 3); % required to find azimuth and elevation of maximums
[~, ~, azim_idx, elev_idx] = find_max_direction(tmp_plane_data, scanning_phi, scanning_theta);
clear tmp_plane_data;
%% find range maximum
[~, r_idx] = max(beams(elev_idx, azim_idx, :));

%% filter maximum
window = 2;

filtered_beams = beams;
filtered_beams(:, :, r_idx - window / 2 : r_idx + window / 2) = 0;

%% Results
all_max = max(abs(beams), [], [3, 2, 1]);

all_filtered_max = max(abs(filtered_beams), [], [3, 2, 1]);
all_filtered_rms = rms(filtered_beams, [3, 2, 1]);

end