% Function to apply HRTFs from Hendrik Kaysers
% HRIR-Database to a given stereo signal
% function [convolved_out] = convolveHRIR(signal,fs,angle)
% Input parameter:
%       signal  = stereo signal
%       fs      = sampling frequency of stereo input signal, should be
%                 48 kHz. This function will resample the impulse response 
%                 to the input sampling frequency,
%                 if the sampling frequency mismatches.
%       angle   = angle of sound incident direction. Currently this will be
%                 limited to -80°:5°:80°, because the HRIRs from anechoic
%                 listening conditions will be applied
% Output parameter:
%       convolved_out = stereo signal with HRIR applied.
%
% Author: Ben Williges
function [convolved_out] = convolveHRIR(signal,fs,angle)
    assert(min(size(signal))>= 2, 'at least two-channel input file needed');
    data = loadHRIR('Anechoic', 80, 0, angle,[,'bte']);
    % resample HRIR (fs = 48kHz) to fs sampling rate (e.g. 44100)
    data.data = resample(data.data,fs,data.fs);
    for ii = 1:size(signal,2) % convolve for each channel
    convolved_out(:,ii) = conv(signal(:,ii),data.data(:,ii));
    end
end