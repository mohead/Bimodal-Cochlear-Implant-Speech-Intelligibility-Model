% Script Test_setSignaltoRMS.m for testing of the function setSignaltoRMS().
clear; close all; clc;
% create random signal a with some rms:
signal_a = randn(512,2);

% create random signal b with some rms:
signal_b = randn(512,2).*3;

% display current rms:
disp('RMS of signal_a:');
rms(signal_a,1)
disp('RMS of signal_b:');
rms(signal_b,1);

% set signal_a to rms of signal_b:
signal_a = setSignaltoRMS(signal_a,rms(signal_b,1),'linear');
disp('RMS of signal_a should be adjusted to signal_b:')
rms(signal_a,1)

% set signal_b to -30 dB FS:
signal_b = setSignaltoRMS(signal_b,-30,'dB');
disp('RMS of signal_b should now be at -30 dB FS:');
rms(signal_b,1)

% test with channels switched:
signal_c = setSignaltoRMS(signal_a',-30,'dB');
disp('RMS of signal_c should now be at -30 dB FS:');
rms(signal_c,2)
