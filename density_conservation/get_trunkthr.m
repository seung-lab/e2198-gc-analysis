function trunk_thr = get_trunkthr(type)
strat_peaks = [];
for i = 1:length(type)
    if ~isempty(str2num(type(i))) && isreal(str2num(type(i)))
        strat_peaks(i) = str2num(type(i));
    end
end
trunk_thr = max(strat_peaks)/10 + 0.1;
end