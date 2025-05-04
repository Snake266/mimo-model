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

%% beamforming
abs_scan_limit = 1/lambda;
scanning_theta = -90 : 1 : 90; 
scanning_phi = -90 : 1 : 90;
    
%% Signal forming
received_signals = zeros(Nrx, size(ps_phases, 2)); % TODO: make this as cell
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
virt_array = cell(Nrx, Ntx);
for rx = 1 : Nrx
    for tx = 1 : Ntx
        corr_val = xcorr(received_signals(rx, :), ps_phases(tx, :), 0);

        virt_array{rx, tx} = corr_val;
    end
end


%% Beamforming
beams = cell(Nrx, Ntx);
for ai = 1:length(scanning_phi)
    for ei = 1:length(scanning_theta)
        azim = scanning_phi(ai);
        elev = scanning_theta(ei);
        [RX, TX] = ndgrid(rx_e, tx_e);
        sv = exp(1j*2*pi/lambda * (RX .* sind(elev) .* cosd(azim) + ...
                                   TX .* sind(elev) .* sind(azim)));
    
        % beams{ai, ei} = sum(received_signal .* sv);
        
        weighted = cellfun(@(v, s) v .* s, virt_array, num2cell(sv), 'UniformOutput', false); % Apply the steeting vector on every element of the array
        beams{ai, ei} = sum(cat(1, weighted{:}), 1); % glue all elements ([Nrx*Ntx, L]) and sum them 
        % filter target maximum
        % max_map(ai, ei) = ...;
        % mean_filt_map(ai, ei) = ... ;
        rp(ai, ei) = mean(abs(beams{ai, ei}).^2); % probably they are reversed (I dunno)
        % rms_filt_map(ai, ei)= ...;
  
    end
end
%% Visualisation
figure;
imagesc(scanning_phi, scanning_theta, 10*log10(rp));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Диаграмма направленности');
colormap('jet');
grid on;

