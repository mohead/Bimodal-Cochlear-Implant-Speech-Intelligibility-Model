function [ output_sig ] = CI_Sim( input_signal,fs,CIType)
%function [ output_sig ] = CI_Sim( input_signal,fs,CIType )
% This function simulates either a Cochlear CI signal processing or a Med El
% CI signal processing. It will firstly generate the electrodogramms using
% scripts, which were sucessfully tested on Med EL subjects, and should work on
% Cochlear ACE subjects.
        def.debug = 0; %Debug-Variable (0 = No Info, 1 = text info, 2 = plotting info)

        CI = createCI();
        CI_Params; %Located in Parameterstore
        %get parameters
        switch CIType
        case 'Cochlear'
           addpath([fileparts(which(mfilename)),filesep,'..',filesep,'CI_CODING_STRAT']); 
           params_CI = CI.setParameter(input_signal,fs,'CI',params_cochlear{:});
           %Ci Cochlear processing left ear
           par=CIParameters_1;
           par.m=length(params_CI.center_frequencies_hz_auralisation);
           par.C=params_CI.MCLl;
           par.T=params_CI.TCLl;
           par.vol=params_CI.Volume;
           TSR = unique(par.pps)*par.n; %Calculate Total stimulation rate = sampling frequency of electrodogramm signal
           par.brand = 0; %Randomization on (1) or off(0)
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % Code part to caliberate the signal to put it in the middle of
           % the dynamic range
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % Set input signal to correct CI calibration level:
           CalibLevel=10; % Base level + x dB [dB] (dynamic = 40 dB)
           DesiredLevel=20*log10(par.B)+CalibLevel;
           desiredAmplitude=10^(DesiredLevel/20);
           inputRMS=rms(input_signal); % Root mean square
           input_signal_CIsim=input_signal/inputRMS*desiredAmplitude;
           % value used for BSIM model:  2.682337941553684e-05
           % Noise from coming the front, with 65 dBFs, using olsa noise file.
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           gain_to_scale = 26.823292988176778;% value used for Linearity test: 3.873933369439492e-4; 
           % values used before caliberation meeting 2.682337941553684e-05 CAUTION: This is the calibration 

           input_signal_CIsim = scaleSignalbyRMS(input_signal,gain_to_scale,'linear');
           [electrodogram, vDistance, StimOrder,vTime, mEnv]=CIStrat(input_signal_CIsim,par,fs);
           params_CI.rms_per_channel = rms2(mEnv,2);
           debug(['Finished ' CIType ' CI processing'],def,@PlotElectrodogram,electrodogram,linspace(0,length(electrodogram)./(TSR),length(electrodogram)),['Electrodogramm ' CIType ' CI']);
           electrodogram_upsampled_ben = resample_elec_zeros(electrodogram,TSR,5*TSR);
           params_CI.voc_sampling_frequency_hz = 5*TSR;
           debug('Finished upsampling of electrodogramm',def,@PlotElectrodogram,electrodogram_upsampled_ben,linspace(0,length(electrodogram_upsampled_ben)./(5*TSR),length(electrodogram_upsampled_ben)),['Electrodogramm ' CIType ' CI upsampled']);
           [output_sig] = CI.Auralisation(electrodogram_upsampled_ben,fs,params_CI);
           %Account for delay
           %output_sig = delay_signal(output_sig,fs,7e-3); %Add delay in miliseconds in front of the signal
        case 'MED-EL'
            params_CI = CI.setParameter(input_signal,fs,'CI',params_med_el{:});
        	[electrodogramm, params_CI] = CI.Simulation(input_signal(:,1),fs,'CI',params_CI);
            params_CI.bandwidth_factor = ones(size(params_CI.bandwidth_factor)); %Set bandwidth for auralization gammatonefilters to 1.
            [output_sig] = CI.Auralisation(electrodogramm,fs,params_CI);
            %Account for delay
            %output_sig = delay_signal(output_sig,fs,7e-3); %Add delay in miliseconds in front of the signal
        otherwise
        error(['This CI Type is not know! Allowed ''Cochlear'', ''Med-El''. You entered: ' CIType]);
        end
        % set output signal to same rms as input signal
        output_sig = setSignaltoRMS(output_sig,rms2(input_signal,1),'linear');
        debug('Finished CI Processing and Auralisation',def,@plotChannels,output_sig',fs,1,['Auralised ' CIType ' CI']);
end

