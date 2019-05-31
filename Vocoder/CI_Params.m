% script for setting vocoder-parameters, which are common across a range of
% experiments.

%Params for Med-el Flex24
params_med_el = {...
'voc_sampling_frequency_hz', 48000,...% Sampling frequency of the electrodogramm
'bandwidth_factor', [1 1 1 1 1 1 1 1 1 1 1 1].*3,... %ERB-Bandwidth of the analysis/synthesis filters. If one electrode needs to be disabled, change the bandwidth accordingly
'gamma_order_stimulation', 3,...
'gamma_order_auralisation', 3,...
'center_frequencies_hz_stimulation', [120 235 384 579 836 1175 1624 2222 3019 4084 5507 7410],...% Middle frequency of analysis-Filterbank
'pulselength', 10e-6,... % Default pulse-length in seconds.
'ipg', 2.1e-6,... % Default inter-pulse-gap in seconds.
'center_frequencies_hz_auralisation', [357 548 689 968 1483 2228 3319 4670 6630 9758 12530 15374],... %[357 548 689 968 1483 2228 3319 4670 6630 9758 12530 15374],... %Based on Data from Landsberger (Insertion angle converted to physical frequency)
'bcompress', 1,... % Switch for compression (0= off, (default), 1 = on)
'B', 0.0156,... %Base level of compression (-34 dB FS) (Laneau Ph.D-Thesis, 2005, p. 124),(Harzcos, Ph.D-Thesis, p. 18).
'M', 1.5859,... %Saturation level of compression (+4 dB FS) (Reference: See above, but for 40 dB input dynamic range (Cochlear Clinical Guidance Document, p. 16)
'Volume', 1,... %Volume control of Med-El-CIs [0 1], to be applied in the CU-conversion step, not directly needed for Vocoder-auralisation
'alpha_c', 340.83,...; %controls the steepness of the compression function. Here we use a value from Stefan Fredelakes ACE implementation. There are other values to be found in the literature, (see Diss Tamas Harzcos, p. 18. for a value of 416.2)
'TCLr', [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]',...
'MCLr', [800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800]',...
'TCLl', [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]',...
'MCLl', [800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800]',...
'spatial_spread', 5,... %After how many electrodes the spatial spread should be neglected (1% of original amplitude remaining at this electrode)
'lambda', 0.0036,... % exponential decay constant for spatial spread (dB/m)
'distance_electrodes', 0.0024,... %Distance between adjunct electrodes [m] for calculation of spatial spread. Estimated from Landberger (2015) Data.
'n_of_m', 12,... %How many electrodes you want to have active at the same time
'pps', 800,... % 800 pps
'electrodeselmethod', 'sequential',...$'random',... % or 'sequential' (fixed order)
'apply_aural_gain', 1,...
'input_gain', 425.5597,... %used for linearity test:  0.006146184930627,... used before caliberatoin for BSIM: 0.000423101334812171,...
'gain_per_channel', 0, ...; %Gain to apply, if apply_aural_gain is set
'debug',0,... %Switches debug mode on (1) or off (0)
'weights',[0.98 0.98 0.98 0.68 0.68 0.45 0.45 0.2 0.2 0.15 0.15 0.15]'}; %weights of analysis-Filterbanks 
  

%Params for Cochlear Contour Advance
params_cochlear = {...
'voc_sampling_frequency_hz', 7200,...% This is the sampling frequency of the ACE-Output (n*maxima * 900 pps). This will be upsampled to 5*7200
'bandwidth_factor', ones(22,1).*1,... %ERB-Bandwidth of the analysis/synthesis filters. If one electrode needs to be disabled, change the bandwidth accordingly
'gamma_order_stimulation', 3,...
'gamma_order_auralisation', 3,...
'center_frequencies_hz_stimulation', [250 375 500 625 750 875 1000 1125 1250 1438 1688 1938 2188 2500 2875 3313 3813 4375 5000 5688 6500 7438],...% Middle frequency of analysis-Filterbank
'pulselength', 10e-6,... % Default pulse-length in seconds.
'ipg', 2.1e-6,... % Default inter-pulse-gap in seconds.
'center_frequencies_hz_auralisation', [722 797 939 1075 1248 1428 1668 1943 2282 2659 3149 3655 4209 4847 5458 6057 7018 8322 9647 11057 12388 13879],...
'bcompress', 1,... % Switch for compression (0= off, (default), 1 = on)
'B', 0.0156,... %Base level of compression (-34 dB FS) (Laneau Ph.D-Thesis, 2005, p. 124),(Harzcos, Ph.D-Thesis, p. 18).
'M', 1.5859,... %Saturation level of compression (+4 dB FS) (Reference: See above, but for 40 dB input dynamic range (Cochlear Clinical Guidance Document, p. 16)
'Volume', 1,... %Volume control of Med-El-CIs [0 1], to be applied in the CU-conversion step, not directly needed for Vocoder-auralisation
'alpha_c', 340.83,...; %controls the steepness of the compression function. Here we use a value from Stefan Fredelakes ACE implementation. There are other values to be found in the literature, (see Diss Tamas Harzcos, p. 18. for a value of 416.2)
'TCLr', [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]',...
'MCLr', [800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800]',...
'TCLl', [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]',...
'MCLl', [800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800]',...
'spatial_spread', 10,... %After how many electrodes the spatial spread should be neglected (1% of original amplitude remaining at this electrode)
'lambda', 0.0036,... % exponential decay constant for spatial spread in dB/m
'distance_electrodes', 0.00099,... %Distance between adjunct electrodes [m] for calculation of spatial spread. Estimated from Landberger (2015) Data.
'n_of_m', 8,... %How many electrodes you want to have active at the same time
'pps', 750,... % 800 pps
'electrodeselmethod', 'sequential',... % or 'sequential' (fixed order)
'apply_aural_gain', 1,...
'gain_per_channel', 0,... %Gain to apply, if apply_aural_gain is set
'debug',0,... %Switches debug mode on (1) or off (0)
'weights',[0.98 0.98 0.98 0.98 0.98 0.98 0.98 0.98 0.98 0.68 0.68 0.68 0.68 0.65 0.65 0.65 0.65 0.65 0.65 0.65 0.65 0.65]'}; %weights of analysis-Filterbanks 

%Params for Acoustic only lowpass
params_acoustic = {...
'voc_sampling_frequency_hz', 48000,...% Bandwidth of Filters
'bandwidth_factor', [1 1 1 1 1 1 1 1 1 1 1 1].*3,... %ERB-Bandwidth of the analysis/synthesis filters. If one electrode needs to be disabled, change the bandwidth accordingly
'gamma_order_stimulation', 3,...
'gamma_order_auralisation', 3,...
'center_frequencies_hz_stimulation', [120 235 384 579 836 1175 1624 2222 3019 4084 5507 7410],...% Middle frequency of analysis-Filterbank
'pulselength', 40e-6,... % Default pulse-length in seconds.
'ipg', 2.1e-6,... % Default inter-pulse-gap in seconds.
'center_frequencies_hz_auralisation', [120 235 384 579 836 1175 1624 2222 3019 4084 5507 7410],...%       
'bcompress', 0,... % Switch for compression (0= off, (default), 1 = on)
'B', 0.0156,... %Base level of compression (-34 dB FS) (Laneau Ph.D-Thesis, 2005, p. 124),(Harzcos, Ph.D-Thesis, p. 18).
'M', 1.5859,... %Saturation level of compression (+4 dB FS) (Reference: See above, but for 40 dB input dynamic range (Cochlear Clinical Guidance Document, p. 16)
'Volume', 1,... %Volume control of Med-El-CIs [0 1], to be applied in the CU-conversion step, not directly needed for Vocoder-auralisation
'alpha_c', 340.83,...; %controls the steepness of the compression function. Here we use a value from Stefan Fredelakes ACE implementation. There are other values to be found in the literature, (see Diss Tamas Harzcos, p. 18. for a value of 416.2)
'TCLr', [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]',...
'MCLr', [800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800]',...
'TCLl', [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]',...
'MCLl', [800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800]',...
'spatial_spread', 5,... %After how many electrodes the spatial spread should be neglected (1% of original amplitude remaining at this electrode)
'lambda', 0.0036,... % exponential decay constant for spatial spread
'distance_electrodes', 0.0028,... %Distance between adjunct electrodes [m] for calculation of spatial spread
'n_of_m', 12,... %How many electrodes you want to have active at the same time
'pps', 800,... % 800 pps
'electrodeselmethod', 'random',... % or 'sequential' (fixed order)
'apply_aural_gain', 0,...
'gain_per_channel', 0,... %Gain to apply, if apply_aural_gain is set
'debug',0,... %Switches debug mode on (1) or off (0)
'weights',[0.98 0.98 0.98 0.68 0.68 0.45 0.45 0.2 0.2 0.15 0.15 0.15]'}; %weights of analysis-Filterbanks 