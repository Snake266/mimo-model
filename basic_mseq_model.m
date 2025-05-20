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
Ntx = 8; Nrx = 8;
tx_e = (0 : (Ntx - 1)) * d;
rx_e = ((0 : (Nrx - 1)) * d)';

%% Generate M-seq
ps_phases = get_mseq_n_times(9, 8);

abs_scan_limit = 1/lambda;
scanning_theta = -25 : 1 : 25;
scanning_phi = -5 : 1 : 45;

beams = mimo_model(ps_phases, tx_e, rx_e, scanning_phi, scanning_theta, targets, lambda);

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



%% Visualisation
figure;
mesh(scanning_phi, scanning_theta, max(abs(beams), [], 3));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Угловые координаты цели (regular max)');
grid on;

figure;
imagesc(scanning_phi, scanning_theta, max(abs(beams), [], 3));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Угловые координаты цели (regular max)');
grid on;

figure;
mesh(scanning_phi, scanning_theta, max(abs(filtered_beams), [], 3));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Отфильтрованные угловые координаты цели (max)');
grid on;

figure;
mesh(scanning_phi, scanning_theta, mean(abs(filtered_beams), 3));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Отфировальные угловые координаты цели (mean)');
grid on;

figure;
mesh(scanning_phi, scanning_theta, rms(abs(filtered_beams), 3));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Отфильтрованные угловые координаты цели (rms)');
grid on;

figure;
mesh(1 : size(filtered_beams, 3), ...
    scanning_theta, ...
    reshape(max(abs(filtered_beams), [], 2), size(filtered_beams, 2), size(filtered_beams, 3)));
xlabel('Дальность (отсчеты)');
ylabel('Угол места (градусы)');
title('Отфильтрованные координаты цели (rms)');
grid on;

%% Results
all_max = max(abs(beams), [], [3, 2, 1])
all_rms = rms(beams, [3, 2, 1])
all_mean = mean(abs(beams), [3, 2, 1])

all_filtered_max = max(abs(filtered_beams), [], [3, 2, 1])
all_filtered_rms = rms(filtered_beams, [3, 2, 1])
all_filtered_mean = mean(abs(filtered_beams), [3, 2, 1])