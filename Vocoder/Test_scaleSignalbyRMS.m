% Script Test_scaleSignalbyRMS.m for testing of the function scaleSignalbyRMS().
clear; close all; clc;
% create random signal a with some rms:
signal_a = randn(512,2);

% display current rms:
disp('RMS of signal_a at beginning:');
rms(signal_a,1)

% scale signal_a by a factor of 2 = + 6 dB:
signal_a = scaleSignalbyRMS(signal_a,2,'linear');
disp('RMS of signal_a should be +6 dB or rms*2:')
rms(signal_a,1)

% scale signal_a by a factor of 0.5 = -6 dB:
signal_a = scaleSignalbyRMS(signal_a,-6,'dB');
disp('RMS of signal_a should now be the same as in the begining:');
rms(signal_a,1)

% test if scaling also works with channels x time signals
signal_c = scaleSignalbyRMS(signal_a',-6,'dB');
disp('RMS of signal_c should now be 6 dB lower again:');
rms(signal_c,2)