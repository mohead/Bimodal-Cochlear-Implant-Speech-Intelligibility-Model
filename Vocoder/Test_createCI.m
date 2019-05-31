% Test script for create ci:
% The test order is from simple to complex.
% Simple functions first, more complex ones (which might call the simple
% functions) later.
% If you have Matlab R2013+ you can run this script using result =runtests('Test_createCI')
% which will make use of the Unittesting framework included since this
% release. If you use an older Matlab version, simply run this script.

clear; close all; clc;

%Test-signal:
[signal,fs] = wavread('OLSA.wav');

% create CI-Object and get all function handles:
CI = createCI;

%% sorting of selected electrode test:
disp('Sequential electrode numbers:');
aa = [12:-1:1]
disp('Now randomized electrode numbers:');
bb = CI.sortElectrodes(aa,'random')
disp('Now they should be in the same order as aa again:');
bb = CI.sortElectrodes(bb,'sequential')
assert(isequal(aa,bb),'Sequential Electrode sorting is faulty. sortElectrodes()');

%% pps-delay calculation test with fixed total pulse rate (12*50 = 600
% pulses per second over all channels, TPR)
hopsize_half = CI.ppshopsize(44100,50,6,12); % Higher pulse rate per electrode = lower hopsize 
hopsize = CI.ppshopsize(44100,50, 12, 12); %fixed total pulse rate over all electrodes, no electrode dropouts = higher hopsize (hopsize_half*2)
assert(hopsize_half == 441,'pps Calculation is faulty');
assert(hopsize == hopsize_half*2,'pps Calculation is faulty');

%% n-of-m electrode selection test:
aa = randn(12,44100*0.5)+1;
hopsize = CI.ppshopsize(44100,10, 4, 12); % TPR = 120 pps
dd = CI.selectElectrodes(aa,4,44100, hopsize);
figure; imagesc(aa); title('Before n of m');
figure; imagesc(dd); title('After n of m'); %Should display only 4 active channels in each 50% overlapping frame

%% Test for correct pps-estimation (ignoring pulse-width):
fs_pps = 6; % ridiculously low
electrodogramm_pps = [1 1 0 1 1 0; ... % 2pps, 2 samples pulse width
                      0 0 1 0 0 1;... % 2 pps, 2 sample pulse width
                      0 1 0 1 0 1;... % 3 pps, 1 sample pulse width [3 pps 3 pps]
                      1 0 1 0 0 1]; % 3 pps, 1 sample pulse width, 2 pps, 1 sample pulse width [3 pps 2 pps]
[pps_hat, ipipps_hat, ipipps_median] = CI.estimate_pps(electrodogramm_pps,fs_pps);
assert(isequal(pps_hat, [2; 2; 3; 3]), 'Error in pps estimation');
assert(isequal(ipipps_hat(~isnan(ipipps_hat)), [2; 2; 3; 3; 3; 2;]), 'Error in interpuls-interval pps estimation'); %Remember that ~isnan makes a VECTOR, thats why the answervektor has this sorting.
assert(isequal(ipipps_median, [2; 2; 3; 2.5]), 'Error in median interpuls pps estimation');

%% Test for perfect reconstruction of n-of-m signal:
aa = randn(12,44100*0.5)+1;
hopsize = CI.ppshopsize(44100,10, 4, 12); % TPR = 120 pps
dd = CI.selectElectrodes(aa,12,44100,hopsize);
bb = [aa  zeros(12,size(dd,2)-size(aa,2))]; %Additional zero-padding for original signal
figure; plot(dd(1,:)./bb(1,:)); %Should equal 1 in every sample
title('Perfect reconstruction = line at 1')
assert(min(dd(1,:)./bb(1,:)) == 1,'Error in perfect reconstruction');
assert(max(dd(1,:)./bb(1,:)) == 1,'Error in perfect reconstruction');

%% assersequentialStimulation-function test:
aa = [eps 0.75 0.25 0; 1 0.25 0.5 1] %numbers test, pulselength 1 sample, keep maxima
bb = [0 1 1 1 0 0 0 0 0;0 0 2 2 2 2 2 2 0;0 1 1 1 0 0 0.5 0.5 0.5;] %shifted puls test, pulselength 3 sample, keep pulses in lowest channel (top one)
dd = CI.testsequentialStimulation(aa);
ee = CI.testsequentialStimulation_fast(bb,3);
expected = [0 0.75 0 0; 1 0 0.5 1];
assert(isequal(dd,expected),'Removal of concurrent stimulation across channels is faulty');
assert(isequal(ee(1,:),[0 1 1 1 0 0 0 0 0]), 'Removal of overlapping channels is faulty');
assert(isequal(ee(2,:),[0 0 0 0 0 2 2 2 0]),'Removal of overlapping channels is faulty');
assert(isequal(min(min(ee(2:3,:))),0),'Removal of overlapping channels is faulty');
% dd should be: dd =
%         0    0.7500         0         0
%    1.0000         0    0.5000    1.0000

%% Test for RIB-2 pulse-timing vector generation
% fake electrodogramm:
elec = zeros(2,10);
elec(1,1:4:end) = 1;
elec(2,2:4:end) = 1;
elec(2,3:4) = 0.8; %One pulse with thrice the pulse duration (and a slightly different amplitude in the "tail", should not increase timing!
[squeezed_electrodogramm, vElectrodes, vDistance] = calculateDistance(elec, 1, 9e-6, 2.1e-6);
assert(isequal(squeezed_electrodogramm, [1 0 1 0 1 0; 0 1 0 1 0 1]), 'Removing zero channels failed!');
assert(isequal(vElectrodes,[1 2 1 2 1 2]),'Wrong electrode selection!');
assert(isequal(round(vDistance*1),[1 1 3 1 3 1]),'Wrong time-distance calculation!');

%% electrodogram-function test:
parameter = CI.setParameter(signal(:,1),fs,'fast_locked','bcompress',1,'debug',1);
% in parameter are now all necessary infos for the CI simulation
disp('fast')
tic
[electrodogramm,parameter] = CI.Simulation(signal(:,1),fs,'fast_locked',parameter);
toc
%zoom in plot
allAxesInFigure = findall(gcf,'type','axes');
set(allAxesInFigure,'XLim',[1.5 1.9])
parameter = CI.setParameter(signal(:,1),fs,'CI','bcompress',1,'electrodeselmethod','random','debug',0);
% in parameter are now all necessary infos for the CI simulation
disp('creation of random sequence');
tic
[electrodogramm_ref,parameter] = CI.Simulation(signal(:,1),fs,'CI',parameter);
toc
parameter = CI.setParameter(signal(:,1),fs,'CI','bcompress',1,'electrodeselmethod','pseudorandom','debug',0);
% in parameter are now all necessary infos for the CI simulation
disp('Applying of random sequence');
tic
[electrodogramm_copy,parameter] = CI.Simulation(signal(:,1),fs,'CI',parameter);
toc
%zoom in plot
allAxesInFigure = findall(gcf,'type','axes');
set(allAxesInFigure,'XLim',[1.5 1.9])

% test for CI-compression and rekompression:
in = [0:0.01:1];
compression = process_compression_ci(in,parameter.B,parameter.M,parameter.alpha_c);
compression_normal = converttoCU(compression,parameter.TCLr(1),parameter.MCLr(1), parameter.Volume);
compression_soft = converttoCU(compression,parameter.TCLr(1),parameter.MCLr(1), parameter.Volume*0.5);
recompressed = inverseCUConversion(compression_normal,parameter.TCLr(1),parameter.MCLr(1), parameter.Volume);
recompressed_soft = inverseCUConversion(compression_soft,parameter.TCLr(1),parameter.MCLr(1), parameter.Volume*0.5);
recompressed = inverseCICompression(recompressed,parameter.B,parameter.M,parameter.alpha_c);
recompressed_soft = inverseCICompression(recompressed_soft,parameter.B,parameter.M,parameter.alpha_c);
figure; 
plot(in,'r');
hold on;
plot(recompressed, 'k');
plot(recompressed_soft, 'k--');
title('orginal signal and reconstructed signal after compression');
assert(max(abs(in(3:59)-recompressed(3:59))) < 5*eps, 'CU compression/recompression is wrong')
assert(max(abs(recompressed-recompressed_soft)) < 5*eps,'Volume conversion in CU compression/recompression is wrong!');

%% Test of gammatone-filter delay:
impulse = zeros(1,44100);
impulse(1,1) = 1;
parameter = CI.setParameter(signal(:,1),fs,'fast_locked','debug',0);
% parameter.center_frequencies_hz_stimulation = [73.2370517960381,107.666366122631,146.019699489555,188.744284893750,236.338327906829,289.356816182397,348.417991091494,414.210556951023,487.501711910993,569.146094147655,660.095747682888,761.411224039868,874.273949189948,1000.00000000000,1140.05545082553,1296.07346920405,1469.87335999932,1663.48178006761,1879.15637082724,2119.41208430771,2387.05050966279,2685.19254212010,3017.31477531529,3387.29004137697,3799.43257149467,4258.54830358123,4769.99092366019,5339.72429446944,5974.39199925237,6681.39481167720;]
% parameter.center_frequencies_hz_auralisation = parameter.center_frequencies_hz_stimulation;
% parameter.voc_sampling_frequency_hz = 16276;
% parameter.bandwidth_factor = 1.*ones(30,1);
% parameter.gamma_order_stimulation = 4;
% parameter.gamma_order_auralisation = 4;
[delayed_impuls_channels, analysis_delay_impulses] = CI.testGammatonedelayline(impulse, parameter, parameter);
figure; 
plotChannels(real(analysis_delay_impulses), 44100,'r');
hold on;
plotChannels(abs(analysis_delay_impulses), 44100, 'r');
title('Output of analysis filter (Envelope and Fine Structure)');
figure; 
plotChannels(delayed_impuls_channels, 44100, 'r');
title('Output of synthethisis channels');

% parameter.center_frequencies_hz_stimulation = [73.2370517960381,107.666366122631,146.019699489555,188.744284893750,236.338327906829,289.356816182397,348.417991091494,414.210556951023,487.501711910993,569.146094147655,660.095747682888,761.411224039868,874.273949189948,1000.00000000000,1140.05545082553,1296.07346920405,1469.87335999932,1663.48178006761,1879.15637082724,2119.41208430771,2387.05050966279,2685.19254212010,3017.31477531529,3387.29004137697,3799.43257149467,4258.54830358123,4769.99092366019,5339.72429446944,5974.39199925237,6681.39481167720;]
% parameter.center_frequencies_hz_auralisation = parameter.center_frequencies_hz_stimulation;
% parameter.voc_sampling_frequency_hz = 16276;
% [delayed_impuls_channels, analysis_delay_impulses] = CI.testGammatonedelayline(impulse, parameter, parameter);

%% Auralization-function test (including delay-line):
parameter = CI.setParameter(signal(:,1),fs,'fast_locked','debug',1);
% in parameter are now all necessary infos for the CI simulation
[electrodogramm,parameter] = CI.Simulation(signal(:,1),fs,'fast_locked',parameter);
vocoded_signal = CI.Auralisation(electrodogramm,fs,parameter);
sound(vocoded_signal,fs);

%% Vocoder-function test:
figure;
vocoded_signal = CI.Vocoder(signal(:,1),fs,'fast_locked', 'debug',1);
sound(vocoded_signal,fs);

% Test different stimulation strategies:
strategies = {'CI', 'Acoustic', 'EAS', 'FSP', 'fast', 'fast_locked'};
for ii = 1:numel(strategies)
    vocoded_signal = CI.Vocoder(signal(:,1),fs,strategies{ii});
end