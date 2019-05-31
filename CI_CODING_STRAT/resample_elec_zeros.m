function [ elec_resampled ] = resample_elec_zeros( elec, old_fs, new_fs)
%function [ elec_resampled ] = resample_elec_ben( elec, old_fs, new_fs)
%   This function resamples an electrodogramm by repeating each sample
%   n-times. Thus it will work only for multiples of the old sampling
%   frequency. It might have the advantage, that an electric pulse is still
%   an electric pulse afterwards. However, there might be additional
%   frequencies in the signal afterwards.
    assert(isequal(size(elec,2),length(elec)),'Input electrodogramm must have a dimension of channels x time-samples');
    n = new_fs/old_fs;
    assert(isequal(n,round(n)),'This function only works for multiples of the original sampling rate');
    elec_resampled = zeros(size(elec,1),size(elec,2)*n);
    elec_resampled(:,1:n:end) = elec; %Just do a n-times zero-padding
end