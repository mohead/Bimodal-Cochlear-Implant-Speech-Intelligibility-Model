% function [valid_pps] = calculate_valid_pps(fs,n,m,pulselength,ipg)
% Function to find valid pps-rates for the CIS strategy given a sampling-rate, pulse-duration
% and n = number of active electrodes out of m (n of m):
% Background on the formula:
% 1) pps*m/n must be an divisor of fs, e.g. 48000 Hz = 60*800 pps *
% (12/12). The divisor is 12. This is for havin an even pps-blocklength in
% samples. The result of this step is the pps_blocklength_samples
% 2) If you have an even pps-blocklength, you must confirm, that all pulses
% can be distributed with equal distance across these blocklength, e.g.
% with at least zero, one or more samples in between your pulses.
% distance_samples = (pps_blocklength_samples - pulselength_samples *n)/n.
% 3) All in one:
% distance_samples = (fs/(pps*n/m) - pulselength_samples *n) / n.
function [valid_pps] = calculate_valid_pps(fs,n,m,pulselength,ipg)
    assert(pulselength < 1e-4,'pulselength must be in seconds and should be in the range of microseconds!');
    assert(ipg < 1e-4, 'inter phase gap must be in seconds and should be in the range of microsenconds!');
    assert(n <= m, 'N must be smaller than or equal to M (N-of-M ACE-like strategy');
    pps = (1:ceil(fs/m)); %All possible pps-rates
    %calculate pulselength:
    pulselength_samples = ceil((2*pulselength+ipg)*fs);
    res = (fs./(pps.*n/m)- pulselength_samples*n)/n;
    res = res(res>=1); % Just use results, which ratio is grater than 1
    valid_pps = pps(res == round(res)); % get all results without any element after decimallimiter.
end