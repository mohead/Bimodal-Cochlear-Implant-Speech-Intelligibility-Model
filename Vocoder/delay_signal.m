function delayed_sig = delay_signal(signal,fs,time)
%This function delays a signal by a given time. It calculates
% the needed samples based on fs and add zeros in front of the signal.
assert(isequal(size(signal,1),length(signal)),'Input signal must be of dimension timesamples x channels!')
assert(round(time*fs)>0,'Delay in seconds must be positive');
samples = round(fs*time);
add_zeros = zeros(samples,min(size(signal)));
delayed_sig = [add_zeros; signal];
end

