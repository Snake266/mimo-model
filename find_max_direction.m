function [est_azim, est_elev, ai_max, ei_max] = find_max_direction(input, azim_array, elev_array)
    [~, max_idx] = max(input(:));
    [ai_max, ei_max] = ind2sub(size(input), max_idx);
    est_azim = azim_array(ai_max);
    est_elev = elev_array(ei_max);
end