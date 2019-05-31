function [ elec_resampled ] = resample_elec_classic( elec, old_fs, new_fs)
%function [ elec_resampled ] = resample_elec_classic( elec, old_fs, new_fs)
%   This function resamples an electrodogramm using the classic upsampling
%   method of adding zeros and then a lowpass-filter for filtering out
%   unwandted aliasing frequencies. This is what is done in matlabs
%   resample.
%   Note that this has the drawback, that pulses are not pulses anymore.
    assert(isequal(size(elec,2),length(elec)),'Input electrodogramm must have a dimension of channels x time-samples');
    elec_resampled = resample(elec',new_fs,old_fs);
    elec_resampled = elec_resampled';
end

