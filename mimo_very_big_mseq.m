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
lambda = 1/fc;
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


ps_phases = reshape(ps, Ntx, []); %% move into function distribute(seq, n_tx);

%% beamforming
abs_scan_limit = 1/lambda;
scanning_theta = -abs : 1 : 90; 
scanning_phi = -90 : 1 : 90;
    
%% Signal forming
received_signals = zeros(Nrx, size(ps_phases, 2));
idx = 1;
for rx = 1 : Nrx
    for tx = 1 : Ntx
        x_tx = tx_e(tx);
        y_rx = rx_e(rx);

        received_phases = exp(1j * 2 * pi / lambda * (x_tx * sind(targets.elevation) * cosd(targets.azimuth) + ...
                                                       y_rx * sind(targets.elevation) * sind(targets.azimuth)));

        received_signals(rx, :) = received_signals(rx, :) + ps_phases(tx, :) * received_phases;
    end
end
%% Virtual array forming
virt_array = cell{Ntx, Nrx};


%% Beamforming
beams = cell{Ntx, Nrx};

for ai = 1:length(scanning_phi)
    for ei = 1:length(scanning_theta)
        azim = scanning_phi(ai);
        elev = scanning_theta(ei);

        sv = exp(1j*2*pi/lambda * (Xv .* sind(elev) .* cosd(azim) + ...
                                   Yv .* sind(elev) .* sind(azim)));
        
        beams{ai, ei} = sum(received_signals .* sv);
        % filter target maximum

        max_map(ai, ei) = ...;
        mean_filt_map(ai, ei) = ...;
        rms_filt_map(ai, ei)= ...;
       
    end
end
power_map = power_map / max(power_map(:));
%% Visualisation
%TODO
figure;
imagesc(scanning_phi, scanning_theta, 10*log10(power_map));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Диаграмма направленности');
colormap('jet');
grid on;
