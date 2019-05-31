classdef CIParameters_1
    
    properties
        m=22; % No. of available electrodes. For Cochlear = 12.
        nDeactivated=[]; % Ch number of deactivated Electrodes. For example electrode 7 == 1 deactivated E.
        n=[8]; % activated electrodes in one cycle. CAVE: If Channels are deactivated, n has to be adjusted: 1 electrode deactivated means n = m-1
        pps=[900*ones(1,22)]; % stimulation rate. Currently only one rate can be chosen for all electrodes
        ipg=2.1*10^-6;; % inter phase gap
        PulseWidth=40; 
        CIFs=16000;
        Sequence=[22:-1:1]; % Order of electrode stimulation
        T
        C
        B=0.0156; % Base Level = level that is mapped to individual T-level 
        M=1.5859; % Maximum level = level that is mapped to individual C-level
        CompFac=340.83; % Compresison Factor that defines steepness of loudness growth function
        vol % volume control limits C-level: 1 = full dynamic range is used
        brand = 0 %randomization of pulse-pattern on (1) or off (0), useful for Auralization part of the vocoder
    end
end
