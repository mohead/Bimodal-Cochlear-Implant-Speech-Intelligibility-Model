function [ output_dB ] = dB_simulations( input,dim )
%[ output_dB ] = dB( input,dim )
% This function returns the sound pressure measured in pascal for the given Input signal.
% It relates an Input of 1 Pascal sound pressure to a sound
% pressure level of 94 dB SPL (exact value: 93.9794 dB SPL), which should correspond to a linear RMS of
% 1. This is also one of the possible reference parameters in the
% AMToolbox, see dbspl.m in there
%     rms(1)= -20*log10(20e-6) = 93.98. This corresponds to the common
%     convention of the reference being 20 micro Pa. Using this
%     convention, the numerical values of signals is the sound pressure
%     measured in Pascal.
output_dB = 120+20*log10((rms2(input,dim))); %eps to prevent from -inf values
end

