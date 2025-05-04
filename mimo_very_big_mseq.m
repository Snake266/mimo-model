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
Ntx = 5; Nrx = 5; % TODO: сейчас количество передатчиков ровно чтобы М 
                % последовательность делилась на числа передатчиков
tx_e = (0 : (Ntx - 1)) * d;
rx_e = ((0 : (Nrx - 1)) * d)';

[Xv, Yv] = meshgrid(tx_e, rx_e');
Xv = Xv(:);
Yv = Yv(:);

%% Generate M-seq
n_mseq = 16;
str = dec2bin(primpoly(n_mseq, 'nodisplay'));
poly = strlength(str)-find(str=='1');
initial = [zeros(1, n_mseq - 2), 1, 1];
pnSequence = comm.PNSequence('Polynomial', poly, ...
    'InitialConditions', initial, ...
    'SamplesPerFrame', 2^n_mseq - 1);
ps = step(pnSequence);
ps = 2 * ps - 1;


ps_phases = reshape(ps, Ntx, []); %% move into function distribute(seq, n_tx);

abs_scan_limit = 1/lambda;
scanning_theta = -45 : 1 : 45; 
scanning_phi = -45 : 1 : 45;
    
%% Signal forming
received_signals = zeros(Nrx, size(ps_phases, 2)); % TODO: make this as cell
for rx = 1 : Nrx
    for tx = 1 : Ntx
        x_tx = tx_e(tx);
        y_rx = rx_e(rx);

        received_phases = exp(1j * 2 * pi / lambda * (x_tx * sind(targets.azimuth) * cosd(targets.elevation) + ...
                                                       y_rx * sind(targets.elevation)));

        received_signals(rx, :) = received_signals(rx, :) + ps_phases(tx, :) * received_phases;
    end
end
%% Virtual array forming
virt_array = cell(Nrx, Ntx);
for rx = 1 : Nrx
    for tx = 1 : Ntx
        corr_val = xcorr(received_signals(rx, :), ps_phases(tx, :), 0);

        virt_array{rx, tx} = corr_val;
    end
end


%% Beamforming
beams = cell(size(scanning_phi, 2), size(scanning_theta, 2));
max_map = zeros(size(scanning_phi, 2), size(scanning_theta, 2));
mean_filt_map = zeros(size(scanning_phi, 2), size(scanning_theta, 2));
rms_filt_map = zeros(size(scanning_phi, 2), size(scanning_theta, 2));
for ai = 1:length(scanning_phi)
    for ei = 1:length(scanning_theta)
        azim = scanning_phi(ai);
        elev = scanning_theta(ei);
        [RX, TX] = ndgrid(rx_e, tx_e);
        sv = exp(-1j*2*pi/lambda * (TX .* sind(azim) .* cosd(elev) + ...
                                   RX .* sind(elev)));
    
        % beams{ai, ei} = sum(received_signal .* sv);
        
        weighted = cellfun(@(v, s) v .* s, virt_array, num2cell(sv), 'UniformOutput', false); % Apply the steeting vector on every element of the array
        beams{ai, ei} = sum(cat(1, weighted{:}), 1); % glue all elements ([Nrx*Ntx, L]) and sum them 
        % filter target maximum
        max_map(ai, ei) = max(abs(beams{ai, ei}).^2);
        mean_filt_map(ai, ei) = mean(abs(beams{ai, ei}.^2));
        rms_filt_map(ai, ei)= rms(abs(beams{ai, ei}.^2));
  
    end
end
%% Visualisation
figure;
imagesc(scanning_phi, scanning_theta, max_map);
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Угловые координаты цели (max)');
grid on;

figure;
imagesc(scanning_phi, scanning_theta, mean_filt_map);
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Угловые координаты цели (mean)');
grid on;


figure;
imagesc(scanning_phi, scanning_theta, rms_filt_map);
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Угловые координаты цели (rms filter)');
grid on;


% virt_amp_vec = zeros(Ntx*Nrx,1);
% for i = 1:Ntx
%     for j = 1:Nrx
%         idx = (j-1)*Ntx + i;
%         virt_amp_vec(idx) = sqrt(sum(virt_array{i, j}.^2)); % можно изменить на sqrt(sum(....^2)) если хотите мощность
%     end
% end

% figure;
% scatter3(Xv, Yv, virt_amp_vec, 50, virt_amp_vec, 'filled');
% xlabel('X (Tx)');
% ylabel('Y (Rx)');
% title('Мощность в виртуальной решетке');
% colorbar;

