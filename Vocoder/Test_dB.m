% Test-script for dB-function
clear; close all; clc;
% create random signal a with some rms:
signal_a = randn(512,2);

% Relate it to dB:
disp('The random Signal has a dB of:')
dB(signal_a,1)

% calculate dB for a linear rms-value of 1 and 0
[dB_out] = dB(1,1);
assert(isequal(dB_out,20*log10(1/20e-6)),'dB to linear RMS-value of 1 is not correct');
[dB_out] = dB(0,1);
assert(isequal(dB_out,20*log10((0/20e-6)+eps)),'dB to linear RMS-value of 0 is not correct');