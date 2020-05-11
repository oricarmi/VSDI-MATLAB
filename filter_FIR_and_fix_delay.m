function [ s_time, s_filteredSignal ] = filter_FIR_and_fix_delay( Nom, Den, t , signal )
% A function that filters the signal and fixes the delay
    filteredSignal = filter(Nom,Den,signal);
    % < -- fix delay from the filtering
    delay = floor(length(Nom)/2); % delay is half of filter order
    s_time = t(1:end-delay); % shifted time vector
    s_filteredSignal = filteredSignal; % shifted ECG signal
    s_filteredSignal(1:delay) = []; % make beginning empty 
    % --->

end

