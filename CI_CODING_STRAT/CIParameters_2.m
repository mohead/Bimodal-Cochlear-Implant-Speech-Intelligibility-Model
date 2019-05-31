classdef CIParameters_2
    
    properties
        m=12; % No. of available electrodes. For Med el = 12.
        nDeactivated=[7 8]; % Ch number of deactivated Electrodes. For example electrode 7 == 1 deactivated E.
        n=[10]; % activated electrodes in one cycle. CAVE: If Channels are deactivated, n has to be adjusted: 1 electrode deactivated means n = m-1
        pps=[500*ones(1,12)]; % stimulation rate. Currently only one rate can be chosen for all electrodes
        Distance; % Wenn eine Rate genutzt wird
        ipg=2.1*10^-6;; % inter phase gap
        PulseWidth=40; 
        CIFs=16000;
        Sequence=[12:-1:1]; % Order of electrode stimulation
        T
        C
        B=0.0156; % Base Level = level that is mapped to individual T-level 
        M=1.5859; % Maximum level = level that is mapped to individual C-level
        CompFac=340.83;; % Compresison Factor that defines steepness of loudness growth function
        vol % volume control limits C-level: 1 = full dynamic range is used
    end
end
