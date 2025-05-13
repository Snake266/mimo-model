function beams = mimo_model(tx_signals, tx_e, rx_e, azims, elevs, targets, lambda)
%
% description.
%
% @since 1.0.0
% @param {type} [name] description.
% @return {cell} beams calculated beams.
% @see dependencies
%

    Ntx = length(tx_e); Nrx = length(rx_e);

    received_signals = zeros(Nrx, size(tx_signals, 2));
    %% Signal forming
    for rx = 1 : Nrx
        for tx = 1 : Ntx
            x_tx = tx_e(tx);
            y_rx = rx_e(rx);

            received_phases = exp(1j * 2 * pi / lambda * (x_tx * sind(targets.azimuth) * cosd(targets.elevation) + ...
                                                            y_rx * sind(targets.elevation)));
            received_signals(rx, :) = received_signals(rx, :) + tx_signals(tx, :) * received_phases;
        end
    end

    %% Virtual array forming
    virt_array = cell(Nrx, Ntx);
    for rx = 1 : Nrx
        for tx = 1 : Ntx
            corr_val = xcorr(received_signals(rx, :), tx_signals(tx, :));
            virt_array{rx, tx} = corr_val;
        end
    end

    %% Beamforming
    beams = cell(length(azims), length(elevs));
    for ai = 1 : length(azims)
        for ei = 1 : length(elevs)
            azim = azims(ai);
            elev = elevs(ei);
            [RX, TX] = ndgrid(rx_e, tx_e);
            sv = exp(-1j*2*pi/lambda * (TX .* sind(azim) .* cosd(elev) + ...
                                        RX .* sind(elev)));
            weighted = cellfun(@(v, s) v .* s, virt_array, num2cell(sv), 'UniformOutput', false);    % apply the steering vector
                                                                                                    % on every element of the array
            
            beams{ai, ei} = sum(cell2mat(weighted)); % squash 8x8 cell into one array
        end
    end
end
