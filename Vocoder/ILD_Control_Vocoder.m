% Test script to exctract ILDs from the noise-files of the vocoder
% before/after:
clear; close all; clc;
% ltfatstart; amtstart;

% [filename filepath] = uigetfile('*.wav', 'Choose an audio file', 'Multiselect', 'on');
% if iscell(filename);  % Sorry, matlab does not return cell type, if
%     else loopend = size(filename,1); end        % it contains only one string...
filename = {'0_nach_Vocoder.wav' '0_vor_Vocoder.wav' '90_nach_Vocoder.wav' '-90_nach_Vocoder.wav' '90_vor_Vocoder.wav' '-90_vor_Vocoder.wav'};
filepath = fullfile(fileparts(which(mfilename)),'..','BSIM',filesep);

%% Get signals from zero deg
[signalref_zero, fs] = audioread(strcat(filepath,filename{2}));
[signalvoc_zero, fs] = audioread(strcat(filepath,filename{1}));

%% Get signals from -90 deg:
[signalref_minin, fs] = audioread(strcat(filepath,filename{6}));
[signalvoc_minin, fs] = audioread(strcat(filepath,filename{4}));


%% Get signals from 90 deg:
[signalref_nin, fs] = audioread(strcat(filepath,filename{5}));
[signalvoc_nin, fs] = audioread(strcat(filepath,filename{3}));

% %% create the reference signals manually:
[signal] = audioread('olnoise.wav');

signalref_minin = convolveHRIR([signal signal],fs,-90);
signalref_zero = convolveHRIR([signal signal],fs,0);
signalref_nin = convolveHRIR([signal signal],fs,90);

% Vocode them:
signalvoc_zero = CI_Sim(signalref_zero(:,2),fs,'MED-EL');
signalvoc_zero = [zeros(size(signalvoc_zero,1),1) signalvoc_zero];
signalvoc_minin = CI_Sim(signalref_minin(:,2),fs,'MED-EL');
signalvoc_minin = [zeros(size(signalvoc_minin,1),1) signalvoc_minin];
signalvoc_nin = CI_Sim(signalref_nin(:,2),fs,'MED-EL');
signalvoc_nin = [zeros(size(signalvoc_nin,1),1) signalvoc_nin];

%% now calculate and plot some ILD-spectra of this signals:
% monaural NH
[Pxx_zero, f] = pwelch(signalref_zero(:,1), 2048, [], 2048, 44100);
[Pxx_minnin, f] = pwelch(signalref_minin(:,1), 2048, [], 2048, 44100);
[Pxx_nin, f] = pwelch(signalref_nin(:,1), 2048, [], 2048, 44100);
% binaural NH
[Pxx_zero_bin_left, f] = pwelch(signalref_zero(:,1), 2048, [], 2048, 44100);
[Pxx_zero_bin_right, f] = pwelch(signalref_zero(:,2), 2048, [], 2048, 44100);
[Pxx_minnin_bin, f] = pwelch(signalref_minin(:,1), 2048, [], 2048, 44100);
[Pxx_nin_bin, f] = pwelch(signalref_nin(:,2), 2048, [], 2048, 44100);

% Plotting of single spectra:
figure;
plot(f,10*log10(Pxx_zero_bin_left),'b'); hold on;
plot(f,10*log10(Pxx_zero_bin_right),'r');
plot(f,10*log10(Pxx_minnin_bin),'b--');
plot(f,10*log10(Pxx_nin_bin),'r--');
legend({'0 deg left','0 deg right','-90 deg binaural left','90 deg binaural right'});
title('Binaural spectra');

figure;
plot(f,10*log10(Pxx_zero),'b'); hold on;
plot(f,10*log10(Pxx_minnin),'b--');
plot(f,10*log10(Pxx_nin),'r--');
legend({'0 deg monaural','-90 deg monaural','90 deg monaural'});
title('Monaural spectra');

% Plotting of ILDs
figure;
plot(f, 10*log10(Pxx_nin)-10*log10(Pxx_zero)); hold on; title('Right ILD');
plot(f, 10*log10(Pxx_nin_bin)-10*log10(Pxx_zero_bin_right),'r');
legend({'Monaural NH','Binaural NH'});
figure;
plot(f, 10*log10(Pxx_minnin)-10*log10(Pxx_zero)); hold on; title('Left ILD');
plot(f, 10*log10(Pxx_minnin_bin)-10*log10(Pxx_zero_bin_left),'r');
legend({'Monaural NH','Binaural NH'});

figure; 
plot(f, 10*log10(Pxx_minnin)-10*log10(Pxx_zero)); hold on; title('monaural ILD');
plot(f, 10*log10(Pxx_nin)-10*log10(Pxx_zero),'r');
legend({'Left ILD','Right ILD'});

figure;
plot(f, 10*log10(Pxx_minnin_bin)-10*log10(Pxx_zero_bin_left)); hold on; title('binaural ILD');
plot(f, 10*log10(Pxx_nin_bin)-10*log10(Pxx_zero_bin_right),'r');
legend({'Left ILD','Right ILD'});

% Now lets do the same with the vocoder processed signals:
% monaural Voc
[Pxxvoc_zero, f] = pwelch(signalvoc_zero(:,2), 2048, [], 2048, 44100);
[Pxxvoc_minnin, f] = pwelch(signalvoc_minin(:,2), 2048, [], 2048, 44100);
[Pxxvoc_nin, f] = pwelch(signalvoc_nin(:,2), 2048, [], 2048, 44100);

% Plotting of spectra:
figure;
plot(f,10*log10(Pxxvoc_zero),'b'); hold on;
plot(f,10*log10(Pxxvoc_minnin),'b--');
plot(f,10*log10(Pxxvoc_nin),'r--');
legend({'0 deg monaural','-90 deg monaural','90 deg monaural'});
title('Monaural spectra Vocoder');

% Plotting of ILDs
figure; 
plot(f, 10*log10(Pxxvoc_minnin)-10*log10(Pxxvoc_zero)); hold on; title('monaural ILD Vocoder');
plot(f, 10*log10(Pxxvoc_nin)-10*log10(Pxxvoc_zero),'r');
legend({'Left ILD','Right ILD'});

% now create binaural signals, e.g. one side with -90 and the other side
% with 0° and calculate ILDs for that signal
% (original signals have one column with silence, e.g. inf ILD)
% sigref_minin = [signalref_minin(:,2) signalref_zero(:,2)];
% sigref_nin = [signalref_zero(:,2) signalref_nin(:,2)];
% 
% sigvoc_minin = [signalvoc_minin(:,2) signalvoc_zero(:,2)];
% sigvoc_nin = [signalvoc_zero(:,2) signalvoc_nin(:,2)];

% %% now extract ILDs and ITDs for each signal:
% % reference;
%  [itd_ref_minin, ild_ref_minin, freqs] = getITDandILD(sigref_minin,fs);
%  [itd_ref_nin, ild_ref_nin, freqs] = getITDandILD(sigref_nin,fs);
% % vocoder:
%  [itd_voc_minin, ild_voc_minin, freqs] = getITDandILD(sigvoc_minin,fs);
%  [itd_voc_nin, ild_voc_nin, freqs] = getITDandILD(sigvoc_nin,fs);
%  
%  % Plot things (reference gets multiplied by -1, due to different ears:
%  figure;
%  plot(freqs,ild_ref_minin,'rx'); hold on; plot(freqs,ild_voc_minin,'bd');
%  title('Left ILD');
%  legend({'Unprocessed','Voc processed'});
%  ylim([-4 14])
%  
%  figure;
%  plot(freqs,ild_ref_nin,'rx'); hold on; plot(freqs,ild_voc_nin,'bd');
%  title('Right ILD');
%  legend({'Unprocessed','Voc processed'});
%  ylim([-4 14])