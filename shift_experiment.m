clear all;
clc; close all;

targets = struct( ...
    'elevation', [15], ...
    'azimuth', [25], ...
    'range', [1000], ...
    'speed', []);

Ntx = 8; Nrx = 8;

scanning_theta = -45 : 1 : 50;
scanning_phi = -34 : 1 : 40;

n = 12;
pack_len = fix((2^n-1)/Ntx);
shift_array = 1 : 10;
L = length(shift_array);
big_max = zeros(1, L);
big_fil_max = zeros(1, L);
big_fil_rms = zeros(1, L);

small_max = zeros(1, L);
small_fil_max = zeros(1, L);
small_fil_rms = zeros(1, L);

parfor shift = 1 : L
    ps_phases = get_mseq_n_times(n, 1);
    ps_phases = circshift(ps_phases, shift);
    new_pss = zeros(Ntx, pack_len);
    for i = 1 : Ntx - 1
        start_index = ((i-1) * pack_len) + 1;
        end_index = start_index + pack_len - 1;
        new_pss(i, :) = ps_phases(start_index : end_index);
    end
    [~, ~, peak, fil_max, fil_rms] = mimo_system(targets, new_pss, Ntx, Nrx, scanning_phi, scanning_theta);
    big_max(shift) = peak;
    big_fil_max(shift) = fil_max;
    big_fil_rms(shift) = fil_rms;

    ps_phases = get_mseq_n_times(n, Ntx);
    ps_phases = circshift(ps_phases, shift);
    [~, ~, peak, fil_max, fil_rms] = mimo_system(targets, new_pss, Ntx, Nrx, scanning_phi, scanning_theta);
    small_max(shift) = peak;
    small_fil_max(shift) = peak / fil_max;
    small_fil_rms(shift) = peak / fil_rms;
end

T = table(big_max', big_fil_max', big_fil_rms', ...
    small_max', small_fil_max', small_fil_rms', ...
    'VariableNames', {'big_max', 'big_fil_max', 'big_fil_rms', 'small_max', 'small_fil_max', 'small_fil_rms'});
writetable(T, 'ShiftExperiment.csv', 'Delimiter', ';', 'QuoteString', 'all')

figure;
hold on; grid on;
scatter3(big_max, big_fil_max, big_fil_rms, 10, 'r', 'DisplayName', 'Big')
scatter3(small_max, small_fil_max, small_fil_rms, 36, 'b', 'DisplayName', 'Small')
legend;
xlabel('MAX')
ylabel('sidelobe')
zlabel('RMS')
view(3);

figure;
hold on; grid on;
scatter3(small_max, small_max / small_fil_max, small_max / small_fil_rms)
xlabel('MAX')
ylabel('sidelobe')
zlabel('RMS')
% figure;
% scatter(big_max, big_fil_rms)
% xlabel('Максимум')
% ylabel('RMS отфильтрованных данных')
% title('Зависимость RMS отфильтрованных данных от максимума')
