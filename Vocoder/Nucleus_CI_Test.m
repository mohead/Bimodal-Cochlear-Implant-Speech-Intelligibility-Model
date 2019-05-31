% Test script for comparison of my CI Signal Processing with 
% Stefan Fredelakes one:
clear all, close all; clc;

% Define center frequencies:
 f_lower = [188 313 438 563 688 813 938 1063 1188 1313 1563 1813 2063 2313 2688 3063 3563 4063 4688 5313 6063 6938];
 f_high = [313 438 563 688 813 938 1063 1188 1313 1563 1813 2063 2313 2688 3063 3563 4063 4688 5313 6063 6938 7938]; 
 
 fc = (f_high-f_lower)/2+f_lower; 
 
 % Get CI function handle
 CI = createCI();
 
 % Signal:
 [signal,fs] = wavread('OLSA.wav');
 
 %Stefans electrodogramm:
 [Stefans_electrodogramm,pps,maxima, ~, ~] = ACE_signal_processing(signal(:,1),fs);
 
 parameter = CI.setParameter(signal(:,1),fs,'CI','center_frequencies_hz_stimulation',fc,...
                                            'MCL', repmat(200,22,1),...
                                            'TCL', repmat(100,22,1),...
                                            'voc_sampling_frequency_hz', 16000,...
                                            'n_of_m', 8,...
                                            'bcompress', 1, ...
                                            'pps', 900,...
                                            'debug',1);
 % Get electrodogramm
 [electrodogramm,parameter] = CI.Simulation(signal(:,1),fs,'CI',parameter);
 %Electrodogramm is currently in CUs, we need it in Iamps to compare with
 % Stefans version:
 electrodogramm = conversion_CL2Iamp(electrodogramm,'nucleus');
 
 % Now compare both electrodogramms:
 figure;
 plotChannels(Stefans_electrodogramm, 'r');
 figure; 
 plotChannels(electrodogramm, 'b');