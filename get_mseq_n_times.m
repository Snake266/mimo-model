function mseqs = get_mseq_n_times(order, n)
    str = dec2bin(primpoly(order, 'all', 'nodisplay'));
    str = str(1:n, :);
    initial = [zeros(1, order - 2), 1, 1];
    mseqs = zeros(n, 2^order - 1);
    for i = 1 : n
        poly = strlength(str(i, :))-find(str(i, :)=='1');
        pnSequence = comm.PNSequence('Polynomial', poly, ...
        'InitialConditions', initial, ...
        'SamplesPerFrame', 2^order - 1);
        ps = step(pnSequence);
        ps = 2 * ps - 1;
        mseqs(i, :) = ps';
    end
end