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
ps_phases = get_mseq_n_times(16, 1);
pack_len = fix((2^16-1)/Ntx);
new_pss = zeros(Ntx, pack_len);
for i = 1 : Ntx - 1
    start_index = ((i-1) * pack_len) + 1;
    end_index = start_index + pack_len - 1;
    new_pss(i, :) = ps_phases(start_index : end_index);
end

abs_scan_limit = 1/lambda;
scanning_theta = -25 : 1 : 25; 
scanning_phi = -25 : 1 : 25;

beams = mimo_model(new_pss, tx_e, rx_e, scanning_phi, scanning_theta, targets, lambda);
data = cellfun(@(x) mean(abs(x .^ 2)), beams, 'UniformOutput', false);
data_mat = cell2mat(data);

[est_azim, est_elev, azim_idx, elev_idx] = find_max_direction(data_mat, scanning_phi, scanning_theta)

window = 100;


% %% Visualisation
figure;
imagesc(scanning_phi, scanning_theta, cell2mat(data));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Угловые координаты цели');
grid on;

