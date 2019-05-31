%changed: input gain
function [stCI] = createCI() 
%% function [stCI] = createCI()
% Simulates Cochlea Implant signal processing using pulsatile sampling of the
% signal envelope in 12 "Electrode channels". It is also possible to simulate
% electroacoustic hearing or acoustic hearing only (via low frequency bandpass).
% After simulation of cochlea implant signal processing the generated sequential
% pulse patterns across channels  is auralized.
%
% In the auralisation stage each pulse gets filtered with a narrow band filter.
% This narrow band filter  can be placed at the aproximated positions of real
% cochlear implant electrode arrays.
%
% Input parameters:
%       None, this function will return function handels, which have their
%       own input parameters-description
% Optional input parameters:
%	All optional input parameters can be evoked with 'key'-'value' pairs.
%   See function setParameters for the details
%
% Output parameters:
%	stCI = struct containing handles to the following functions:
%     stCI.Vocoder = full CI simulation and auralization for one channel
%     stCI.Vocoder_stereo = full stereo CI simulation and auralization
%     stCI.Simulation = only CI simulation, will return an electrodogramm
%     stCI.Auralisation = takes an electrodogramm as input and will
%     auralize it. You will have to know the bandwise rms and overall rms afterwards.
%     stCI.setParameter = Define Parameters for the CI Simulation.
%     stCI.selectElectrodes = for testing of the n-of-m Electrode selection
%     stCI.testsequentialStimulation = for testing of the
%     sequentialStimulation-function
%
% Author : Ben Williges <ben.williges@uni-oldenburg.de>
%%=======================================================================%%

% add necessary folders to MATLABs PATH
addpath([fileparts(which(mfilename)),filesep,'Gammatonefilterbank']);

stCI.Vocoder = @CI_Vocoder;
stCI.Simulation = @CI_Simulation;
stCI.Auralisation = @CI_auralisation;
stCI.setParameter = @setParameter;
stCI.selectElectrodes = @choose_n_of_m_channels;
stCI.sortElectrodes = @sortElectrodesperBlock;
stCI.ppshopsize = @calculate_delay;
stCI.testsequentialStimulation = @assertsequentialStimulation;
stCI.testsequentialStimulation_fast = @assertsequentialStimulation_fast;
stCI.testGammatonedelayline = @test_gammatonefilterbankdelay;
stCI.testGammatonedelayline_rms = @test_gammatonefilterbankdelay_rms;
stCI.estimate_pps = @estimate_ppsinterval;
end
    function [vocoded_signal, parameter] = CI_Vocoder(signal,fs,vocoder_type, varargin)
        parameter = setParameter(signal,fs,vocoder_type,varargin{:}); %Validates input;
        %% Create Electric stimulation pattern
        [electrodogramm, parameter] = CI_Simulation(signal,fs, vocoder_type ,parameter);
        
        % Add vocoder-type to parameter struct (compability issue)
        parameter.vocodertype = vocoder_type;
        %% CI Auralisation
        [vocoded_signal, parameter] =  CI_auralisation(electrodogramm,fs,parameter);
    end

%% Begin of main functions
    function [electrodogramm, CISIM_parameter] = CI_Simulation(signal,fs, vocoder_type, CISIM_parameter)
        w     = 2*1200/fs;
        [b,a] = butter(1,w,'high');
        signal     = filter(b,a,signal); %Preemphasis filter, see Diss Laneau, 2006, p. 51
        % Resampling to vocoder-frequency
        if fs ~= CISIM_parameter.voc_sampling_frequency_hz;
            signal = resample(signal,CISIM_parameter.voc_sampling_frequency_hz,fs);
        end
        
        % Create analyse -Filterbank
        analyzer_stim = Gfb_Analyzer_new(CISIM_parameter.voc_sampling_frequency_hz,...
            CISIM_parameter.gamma_order_stimulation, ...
            CISIM_parameter.center_frequencies_hz_stimulation,...
            CISIM_parameter.bandwidth_factor);
        
        % First analysis-filterbank
        %1st analyse:
        [y1, analyzer_stim] = Gfb_Analyzer_process(analyzer_stim,signal');
        
        %CI simulation:
        
        % Calculate envelope after first filterbank
        env_signal_new = abs(y1); %Envelope
        weights = repmat(CISIM_parameter.weights,1,size(env_signal_new,2));
        env_signal_new = sqrt(weights.*(real(y1).^2+imag(y1).^2));
        fine_signal_new = real(y1); %Fine structure information
        
        CISIM_parameter.env = env_signal_new;
                
        % calculate rms per analysis-channel:
        CISIM_parameter.rms_per_channel = rms2(y1,2);
        
        % calculate pps (not all strategies might need that)
        CISIM_parameter.block_delay = calculate_delay(CISIM_parameter.voc_sampling_frequency_hz,CISIM_parameter.pps, CISIM_parameter.n_of_m, length(CISIM_parameter.center_frequencies_hz_stimulation));
        %calculate pulse-length in samples:
        CISIM_parameter.len_pulse = ceil((2*CISIM_parameter.pulselength+CISIM_parameter.ipg)*CISIM_parameter.voc_sampling_frequency_hz); %Round towards infinity to make sure, that exact timing is always possible
        
        % Apply n- of m-processing of data_channels if needed
        env_signal_new = choose_n_of_m_channels(env_signal_new, CISIM_parameter.n_of_m, CISIM_parameter.voc_sampling_frequency_hz, CISIM_parameter.block_delay);
        fine_signal_new = choose_n_of_m_channels(fine_signal_new, CISIM_parameter.n_of_m, CISIM_parameter.voc_sampling_frequency_hz, CISIM_parameter.block_delay);
        
        % lowpass filter envelope:
        env_signal_new = lowpass_filter_env(env_signal_new, CISIM_parameter.center_frequencies_hz_stimulation, CISIM_parameter.voc_sampling_frequency_hz);
        CISIM_parameter.env_lp = env_signal_new;
        
        % sample envelope or fine structure depending on coding strategy
        [electrodogramm, CISIM_parameter] = coding_strategy(env_signal_new, fine_signal_new, vocoder_type,CISIM_parameter);
        
        % determine global rms of envelope and fit that in the middle of
        % the dynamic range:
        if CISIM_parameter.apply_aural_gain
           electrodogramm = scaleSignalbyRMS(electrodogramm, CISIM_parameter.input_gain,'linear');
        else
            rms_env_global = rms2(electrodogramm(electrodogramm>0),1); %Note env_signal_new>0 creates a vector over all channels!
            rms_goal = 0.05.*(CISIM_parameter.M-CISIM_parameter.B)+CISIM_parameter.B;
            CISIM_parameter.input_gain = rms_goal/rms_env_global;
            electrodogramm = scaleSignalbyRMS(electrodogramm, CISIM_parameter.input_gain,'linear');
        end
        if CISIM_parameter.debug
            timeAxis    = (1:length(electrodogramm))./CISIM_parameter.voc_sampling_frequency_hz;
            figure(1)
            subplot(1,2,1)
            title('Electrogram Pre vs Post compression')            
            plot (timeAxis,electrodogramm(1,:))
            xlabel('Time (s)')
            ylabel('Amplitude')
        end
        % Compression of pulses and mapping to Clinical Units
        if CISIM_parameter.bcompress %Compress envelope if needed
            electrodogramm = process_compression_ci(electrodogramm,CISIM_parameter.B,CISIM_parameter.M,CISIM_parameter.alpha_c);
            electrodogramm = converttoCU(electrodogramm, CISIM_parameter.TCLr, CISIM_parameter.MCLr, CISIM_parameter.Volume);
        end
        if CISIM_parameter.debug
            subplot(1,2,2)
            plot (timeAxis,electrodogramm(1,:),'r')
            xlabel('Time (s)')
            ylabel('Amplitude in clincal current units')
        end
        % pps-estimation
        CISIM_parameter.pps_channels = estimate_ppsinterval(electrodogramm, CISIM_parameter.voc_sampling_frequency_hz);
        %debug: plot electrodogramm:
        if CISIM_parameter.debug
            electrodogramm_plot = electrodogramm;
            % fine structure
            figure(2)
            plotChannels(fine_signal_new, CISIM_parameter.voc_sampling_frequency_hz, 1,'Finestructure');
            xlabel('Time [s]');
            ylabel('Channel');
            ylim([0.8 12.5]);
            % original envelope signal
            plotChannels(env_signal_new, CISIM_parameter.voc_sampling_frequency_hz, 1, 'Envelope');
            xlabel('Time [s]');
            ylabel('Channel');
            % sampled data
            if CISIM_parameter.bcompress 
                            plotChannels(electrodogramm_plot, CISIM_parameter.voc_sampling_frequency_hz, 1/1200, strcat('Coding strategy: ', vocoder_type));
                            xlabel('Time [s]');
                            ylabel('Elektrode');
                            ylim([0.8 12.5]);
            else
                plotChannels(electrodogramm_plot, CISIM_parameter.voc_sampling_frequency_hz, 1, strcat('Coding strategy: ', vocoder_type));
                xlabel('Time [s]');
                ylabel('Elektrode');
                ylim([0.8 12.5]);
            end
        end
        
    end

    function [outsig, auralisation_parameter] = CI_auralisation(electrodogramm,fs,auralisation_parameter)
        % convert CUs back to normal amplitudes:
        if auralisation_parameter.bcompress
            electrodogramm = inverseCUConversion(electrodogramm,auralisation_parameter.TCLr, ...
                auralisation_parameter.MCLr, auralisation_parameter.Volume);
            electrodogramm = inverseCICompression(electrodogramm,auralisation_parameter.B, ...
                auralisation_parameter.M, ...
                auralisation_parameter.alpha_c);
        end
        % Simulate spatial spread
        % Interaction between electrodes
        if ~any(strcmp(auralisation_parameter.vocoder_type, {'EAS','Acoustic'}))
            [electrodogramm,mspatial_weights] = CIinteract_Ben(electrodogramm, auralisation_parameter.spatial_spread, auralisation_parameter.distance_electrodes, auralisation_parameter.lambda);
        elseif strcmp(auralisation_parameter.vocoder_type,'EAS')%No spatial spread in the acoustic low frequency part
            [electrodogramm(4:end,:),mspatial_weights] = CIinteract_Ben(electrodogramm(4:end,:), auralisation_parameter.spatial_spread, auralisation_parameter.distance_electrodes, auralisation_parameter.lambda);
        end
        % Create analyse -Filterbank
        analyzer_aural = Gfb_Analyzer_new(auralisation_parameter.voc_sampling_frequency_hz,...
            auralisation_parameter.gamma_order_auralisation,...
            auralisation_parameter.center_frequencies_hz_auralisation,...
            auralisation_parameter.bandwidth_factor);
        % Next Analyse-Filterbank (modulates puls-pattern to desired frequency region on basiliar
        % membrane)
        electrodogramm = Gfb_Analyzer_process(analyzer_aural, electrodogramm);
        % Scale to correct channel-rms
        if auralisation_parameter.apply_aural_gain
            electrodogramm = scaleSignalbyRMS(electrodogramm, auralisation_parameter.rms_per_channel, 'linear');
        else
            [electrodogramm, auralisation_parameter.gain_per_channel] = setSignaltoRMS(electrodogramm,auralisation_parameter.rms_per_channel,'linear');
        end
        % Signal-Synthesis:
        synthesizer_aural = Gfb_Synthesizer_new (analyzer_aural, 1/100); %1/100 = Delay
        if auralisation_parameter.debug
            delay_output = Gfb_Delay_process(synthesizer_aural.delay, electrodogramm);
            plotChannels(delay_output,auralisation_parameter.voc_sampling_frequency_hz,1, 'r');
            title('After synthesis delay');
            % calculate and plot frequency response of auralisation filters:
            impulse = zeros(1*auralisation_parameter.voc_sampling_frequency_hz,1); impulse(1) = 1; %creates direc impulse
            % Create analyse -Filterbank for auralisation
            analyzer_aural = Gfb_Analyzer_new(auralisation_parameter.voc_sampling_frequency_hz,...
                                        auralisation_parameter.gamma_order_auralisation,...
                                        auralisation_parameter.center_frequencies_hz_auralisation,...
                                        auralisation_parameter.bandwidth_factor);
            % Next Analyse-Filterbank (modulates puls-pattern to desired frequency region on basiliar 
            % membrane)
            impulse_out = Gfb_Analyzer_process(analyzer_aural, impulse);
            frequency_response = fft(real(impulse_out)');                     
            frequency = [0:length(impulse)-1] * auralisation_parameter.voc_sampling_frequency_hz / length(impulse);                        
            Gfb_plot(160, [0, auralisation_parameter.voc_sampling_frequency_hz/2, -40, 0], ...
            strcat('frequency response of the auralisation filterbank with spatial spread constant lambda [m] = ', num2str(auralisation_parameter.lambda)), ...
            'frequency / Hz', 'filter response / dB', ...
            frequency, 20 * log10(abs(frequency_response)));
            hold on; plot((1:length(frequency)),-3.*ones(size(frequency)),'k-');
            disp(' ');
            disp('Figure 160 shows the frequency response of the auralisation filterbank.');
            % Also add the spatial spread to the figure:
            plot(repmat(auralisation_parameter.center_frequencies_hz_auralisation,size(electrodogramm,1),1)',20*log10(mspatial_weights)','x--','Markersize',15);
        end
        outsig = Gfb_Synthesizer_process(synthesizer_aural, electrodogramm);
        outsig = outsig'; %Make mono channel
        % Resample to original sampling frequency
        outsig = resample(outsig,fs,auralisation_parameter.voc_sampling_frequency_hz);
        val = max(abs(outsig));
        if val>1,
            warning('Signal clipped!');
        end
    end

    function parameter = setParameter(signal,fs,vocoder_type,varargin)
        % Input parser (make all parameters except the first three optional)
        p = inputParser; %creates Input parser structure
        p.addRequired('signal', @(x) isnumeric(x));
        p.addRequired('fs', @isnumeric);
        p.addRequired('vocoder_type', @ischar); % Can be any acceptable switch-case from coding_strategy()
        % PPS-Rate over all channels. To obtain PPS per Channel divide this number by
        % number of center_frequencies_hz_stimulation
        p.addParamValue('voc_sampling_frequency_hz', 48000, @isnumeric);
        % Bandwidth of Filters
        p.addParamValue('bandwidth_factor', 3.*[1 1 1 1 1 1 1 1 1 1 1 1], @isnumeric); %ERB-Bandwidth of the analysis/synthesis filters. If one electrode needs to be disabled, change the bandwidth accordingly
        % Order of Gammatonefilterbank (stimulation and auralisation analysis
        % filterbank)
        p.addParamValue('gamma_order_stimulation', 3, @isnumeric);
        p.addParamValue('gamma_order_auralisation', 3, @isnumeric);
        % Middle frequency of analysis-Filterbank
        p.addParamValue('center_frequencies_hz_stimulation', [120 235 384 579 836 1175 1624 2222 3019 4084 5507 7410], @isnumeric);
        p.addParamValue('pulselength', 40e-6, @isnumeric); % Default pulse-length in seconds.
        p.addParamValue('ipg', 2.1e-6, @isnumeric); % Default inter-pulse-gap in seconds.
        % Higher version with 12 Electrodes using greenwoods equation and an
        % electrode insertion depth of 24 mm in a 32 mm cochlear duct (Flex24 Electrode)
        % old auralization frequencies *with* included dominance region
        % between 500 Hz and 1.2 kHz [390 550 759 1031 1384 1843 2440 3216 4225 5537 7243 9460]
        p.addParamValue('center_frequencies_hz_auralisation', [120 235 384 1794 2181 2646 3203 3872 4675 5638 6793 8179], @isnumeric);
        p.addParamValue('bcompress', 0, @isnumeric); % Switch for compression (0= off (default), 1 = on)
        p.addParamValue('B', 0.0156, @isnumeric); %Base level of compression (see Diss Tamas Harzcos, p. 18 (4/256)
        p.addParamValue('M', 0.5859, @isnumeric); %Saturation level of compression (see Diss Tamas Harzcos, p. 18 (150/256)
        p.addParamValue('Volume', 1, @isnumeric); %Volume control of Med-El-CIs [0 1], to be applied in the CU-conversion step, not directly needed for Vocoder-auralisation
        p.addParamValue('alpha_c', 416.2, @isnumeric); %controls the steepness of the compression function (see Diss Tamas Harzcos, p. 18
        p.addParamValue('TCLr', [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]', @isnumeric);
        p.addParamValue('MCLr', [800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800]', @isnumeric);
        p.addParamValue('TCLl', [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]', @isnumeric);
        p.addParamValue('MCLl', [800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800]', @isnumeric);
        p.addParamValue('spatial_spread', 5, @isnumeric); %After how many electrodes the spatial spread should be neglected (1% of original amplitude remaining at this electrode)
        p.addParamValue('lambda', 0.0036, @isnumeric); % exponential decay constant for spatial spread
        p.addParamValue('distance_electrodes', 0.0028, @isnumeric); %Distance between adjunct electrodes [m] for calculation of spatial spread
        p.addParamValue('n_of_m', 12, @isnumeric); %How many electrodes you want to have active at the same time
        p.addParamValue('pps', 800, @isnumeric);
        p.addParamValue('electrodeselmethod', 'random', @ischar); % or 'sequential' (fixed order)
        p.addParamValue('apply_aural_gain', 0, @isnumeric);
        p.addParamValue('gain_per_channel', 0, @isnumeric); %Gain to apply per vocoder-auralisation-channel, if apply_aural_gain is set
        p.addParamValue('input_gain', 0, @isnumeric);
        p.addParamValue('debug',0, @isnumeric);
        p.addParamValue('weights',[1 1 1 1 1 1 1 1 1 1 1 1]', @isnumeric)
        % Validation
        p.parse(signal, fs, vocoder_type, varargin{:});
        parameter = p.Results;
        % One additional validation:
        if parameter.n_of_m >size(parameter.center_frequencies_hz_stimulation,2)
            error('n_of_m must be equal or smaller than the center frequencies of the CI stimulation (''center_frequencies_hz_stimulation''-Parameter)');
        end
        assert(isequal(length(parameter.bandwidth_factor),length(parameter.center_frequencies_hz_stimulation)),'You must specify an Gammatone ERB-bandwidth for each analyzing channel');
    end

    function [electrodogramm, parameter] = coding_strategy(envelope, finestructure, vocoder_type,parameter)
        switch vocoder_type
            case 'CI'
                electrodogramm = pulsatile_sampling_CIS(envelope,parameter.len_pulse, parameter.block_delay, parameter.electrodeselmethod,parameter.n_of_m);
                electrodogramm = assertsequentialStimulation(electrodogramm);
            case 'EAS'
                electrodogramm = zeros(size(envelope));
                startkanal = 4; % Using CIs-pulspattern from this channel onwards
                electrodogramm(startkanal:end,:) = pulsatile_sampling_CIS(envelope(startkanal:end,:),parameter.len_pulse,parameter.block_delay, parameter.electrodeselmethod,parameter.n_of_m);
                electrodogramm(startkanal:end,:) = assertsequentialStimulation(electrodogramm(startkanal:end,:));
                electrodogramm(1:startkanal-1,:) = finestructure(1:startkanal-1,:); %Replace envelope with real-Part of original signal
                parameter.center_frequencies_hz_auralisation(1:startkanal-1) = parameter.center_frequencies_hz_stimulation(1:startkanal-1); % First auralisation bandpass must be the same as in stimulation
            case 'Acoustic'
                electrodogramm = zeros(size(envelope));
                startkanal = 3;
                % Only use the first three channels
                electrodogramm(1:startkanal,:) = finestructure(1:startkanal,:);
                parameter.rms_per_channel(startkanal+1:end) = 0; %Set rms to zero for all other frequency channels
                parameter.center_frequencies_hz_auralisation(1:startkanal) = parameter.center_frequencies_hz_stimulation(1:startkanal); % First auralisation bandpass must be the same as in stimulation
            case 'FSP'
                electrodogramm = pulsatile_sampling_CIS(finestructure, parameter.len_pulse, parameter.block_delay, parameter.electrodeselmethod,parameter.n_of_m);
                electrodogramm = assertsequentialStimulation(electrodogramm);
            case 'fast'
                electrodogramm = fast(envelope,parameter.len_pulse);
                electrodogramm = assertsequentialStimulation_fast(electrodogramm,parameter.len_pulse);
            case 'fast_locked'
                electrodogramm = zeros(size(envelope));
                startkanal = 4; % Using FAST-Envelope from this channel onwards
                electrodogramm(startkanal:end,:) = fast(envelope(startkanal:end,:),parameter.len_pulse);
                electrodogramm(1:startkanal-1,:) = fast(finestructure(1:startkanal-1,:),parameter.len_pulse);
                electrodogramm = assertsequentialStimulation_fast(electrodogramm,parameter.len_pulse);
            otherwise
                error('This Vocoder configuration is not known');
        end
    end

    function [sampled_puls_pattern, rand_channels] = pulsatile_sampling_CIS(signal_to_sample,pulselength_samples, pps_blockshift, sortmethod, n_of_m)
        % Create pulsatile pattern with correct sampling pulse and randomized sequence between electrodes
        sampled_puls_pattern = zeros(size(signal_to_sample));
        pulse = ones(1,pulselength_samples); % Generate pulse. Note that this pulselength contains negative Phase, inter-phase-gap and positive phase!
        n_samples = 1;
        min_value = 0; 
        nrchans = size(signal_to_sample,1);
        if length(pps_blockshift) > 1
            for ii = 1:length(pps_blockshift)
                channelshift(ii) = round((pps_blockshift(ii)-pulselength_samples*n_of_m)/n_of_m);
            end
        else
        channelshift = round((pps_blockshift-pulselength_samples*n_of_m)/n_of_m);    
        end
        chanshift_idx = 1;
        pps_idx = 1;
        if n_of_m*pulselength_samples > min(pps_blockshift)
            error('Using this pps, N-of-M, pulselength, and fs would result in parallel stimulation! Please adjust these three factors, so that the followind equation is satisfied, with x = [1,2,3,4,...,100]: ((x * N) + pulselength_samples)*pps*M/N = fs');
        end
        if strcmp(sortmethod,'pseudorandom')
            load('Rand_channels.mat');
            
        else
        rand_channels = [];
        end
        rand_chan_idx = 1;
        size_signal = size(signal_to_sample,2);
        while n_samples <= size_signal-1;
                    if strcmp(sortmethod,'pseudorandom')
                        rand_channel = rand_channels(rand_chan_idx,:);
                    else
                        rand_channel = sortElectrodesperBlock([1:size(signal_to_sample,1)], sortmethod);
                    end
                    for n_kanaele = 1:nrchans;
                        %Break before last sample 
                        if n_samples > length(signal_to_sample)-pulselength_samples;
                            break
                        end
                        % Pulsatile sampling if amplitude of finestructure > Threshold (currently 0)
                        if signal_to_sample(rand_channel(n_kanaele), n_samples) > min_value;
                            sampled_puls_pattern(rand_channel(n_kanaele),[n_samples : n_samples+pulselength_samples-1]) = ...
                                signal_to_sample(rand_channel(n_kanaele),n_samples).*pulse; %sampling with puls
                        end
                        if signal_to_sample(rand_channel(n_kanaele), n_samples) > min_value && n_kanaele < size(signal_to_sample,1)
                            n_samples = n_samples+pulselength_samples+channelshift(chanshift_idx); % Shift to next channel (sequential stimulation)
                            if chanshift_idx == length(chanshift_idx)
                                chanshift_idx = 0;
                            end
                            chanshift_idx = chanshift_idx + 1;
                        end
                        if n_samples == length(signal_to_sample)-1;
                            n_samples = n_samples+1;
                            break
                        end
                    end
                    n_samples = n_samples+pps_blockshift(pps_idx)-(n_of_m-1)*(pulselength_samples+channelshift(chanshift_idx)); %Forwards to next block (according to pps)
                    if pps_idx == length(pps_blockshift)
                        pps_idx = 0;
                    end
                    pps_idx = pps_idx +1;
                    rand_chan_idx = rand_chan_idx+1;
                    rand_channels = [rand_channels; rand_channel];
        end
        
        if strcmp(sortmethod,'random')
            save('Rand_channels.mat','rand_channels');
        end
end
    function env_signal_filtered = lowpass_filter_env(env_signal_new, center_frequencies,fs)
        env_signal_filtered = nan(size(env_signal_new));
        for ii = 1:size(env_signal_filtered,1)
            % design 1rd order low-pass butterworth IIR-Filter
            % Envelope filter coeffs with 200 Hz cutoff-frequency
            [b_env, a_env] = butter(1, (200./(fs/2)));
            env_signal_filtered(ii,:) = filtfilt(b_env, a_env, env_signal_new(ii,:));
        end
    end

    function sparse_electrode_pattern = choose_n_of_m_channels(channelgramm, n_of_m, fs, hopsize)
        if size(channelgramm,1) ~= n_of_m %Only do this time-consuming task, if we have to do it
            % Begin block processing
            % Split audio signal into blocks
            windowname = 'rectwin'; % Use rectangular Window
            length_block=min(hopsize); %Number of samples per block
            if rem(length_block,2) ~=0 %check if length block is uneven number and make it even
                length_block = length_block+1;
            end
            window_syn_analyse = window(windowname, length_block);
            
            reconstructed_signal = zeros(size(channelgramm));
            
            startindex=1;
            anz_bloecke=floor((length(channelgramm)-length_block)/length_block)+1; %Number of Blocks
            block_index=[startindex length_block];
            
            for ii=1:anz_bloecke                   %Begin block processing
                % Get current block
                samples=channelgramm(:,block_index(1):block_index(2));
                % Windowing:
                samples = bsxfun(@times,samples,sqrt(window_syn_analyse)');
                %get current maximum channels:
                [~,idx_max_channels] = sort(rms2(samples,2),'descend');
                % Only keep maximum channels:
                samples(idx_max_channels(n_of_m+1:end),:) = 0;
                
                %% Re-synthesis of block processing
                reconstructed_signal(:,block_index(1):block_index(2)) = reconstructed_signal(:,block_index(1):block_index(2)) + (bsxfun(@times,samples,sqrt(window_syn_analyse)')); %Overlap add
                % move one block forwards:
                block_index = block_index+(length_block); % New Blockindex (taking in account the overlap shift)
                if block_index(2) > size(channelgramm,2) % Zero-Padding in last block
                    add_zeros = block_index(2)-size(channelgramm,2);
                    add_zeros = zeros(size(channelgramm,1),add_zeros);
                    reconstructed_signal = [reconstructed_signal add_zeros];
                    break; % This break is important! Without it, we will enter the loop one last time, which will screw up the block processing
                end
            end
            % throw out additional zeroes, they migh screw up the resample of
            % ITDs (happens in one of 120 OLSA-Sentences, unfortunately...)
            reconstructed_signal = reconstructed_signal(1:size(channelgramm,1),1:size(channelgramm,2));
            sparse_electrode_pattern = reconstructed_signal;
        else
            sparse_electrode_pattern = channelgramm;
        end
    end

    function [melectrodenumbers] = sortElectrodesperBlock(electrodenrs,choosemethod)
        switch choosemethod
            case 'random'
                rand_idx = randperm(length(electrodenrs));
                melectrodenumbers = electrodenrs(rand_idx);
            case 'sequential'
                melectrodenumbers = sort(electrodenrs,'descend');
            otherwise
                error('You must select either ''random'' or ''sequential'' for ordering the electrodes per block!');
        end
    end

    function sampled_puls_pattern = fast(signal,len_pulse)
        % General Idea: Find all indexes per Channel, where the envelope has maxima.
        % Place a pulse at the with an amplitude of the maximum
        % at that place.
        sampled_puls_pattern = zeros(size(signal));
        pulse = ones(1,len_pulse); % Generate pulse
        
        for ii = 1:size(signal,1) %For all channels
            % Gets indexes of maxima, do not include endpoints
            sig = 0; % most sensitive peak picking possible
            thres = 1e-10; % most sensitive peak picking possible
            [peakLoc, peakMag] = peakfinder(signal(ii,:),sig,thres,1,0);
            diff_peak = diff(peakLoc); %Get difference in samples betweeen peaks
            %Throw out all peaks which are less than length(pulse) away
            peakLoc = peakLoc(diff_peak>=length(pulse));
            % Apply maximum amplitude to electric puls and place it at the
            % positive slope
            for kk = 1:length(peakLoc)
                sampled_puls_pattern(ii,[peakLoc(kk):peakLoc(kk)+len_pulse-1]) = signal(ii,peakLoc(kk)).*pulse;
            end
        end
        
    end

    function electrodogramm = assertsequentialStimulation(electrodogramm)
        % create temporal matrix with just zeros (no stimulus) or one (stimulus)
        stimpattern = zeros(size(electrodogramm));
        stimpattern(electrodogramm>0) = 1;
        % Calculate sum:
        nr_chans_stimulated = sum(stimpattern,1);
        if any(nr_chans_stimulated>1)
            warning('Parallel electrode stimulation detected! Will only keep the electrode with maximum amplitude!');
        end
        [~,tbins_stimulation] = find(nr_chans_stimulated>1);
        for dd = 1:length(electrodogramm(nr_chans_stimulated>1))
            slice = electrodogramm(:,tbins_stimulation(dd));
            [max_val, max_idx] = max(slice);
            slice(:) = 0; %Set all channels which are stimulated at the same time to zero
            slice(max_idx) = max_val; %keep only the maximum one
            electrodogramm(:,tbins_stimulation(dd)) = slice;
        end
    end

    function electrodogramm_out = assertsequentialStimulation_fast(electrodogramm,pulselen)
        % Assert sequential stimulation for the fast coding strategy
        % Keep pulses, when there are none in the lowest channels,
        % otherwise prefer the lowest channels and discard higher frequency
        % channels
        stimpattern = electrodogramm;
        stimpattern(stimpattern>0) = 1; %Binary mask
        blacklist_idx = find(stimpattern(1,:)>0); %Blacklist of first channel
        for ii = 2:size(stimpattern,1)
            chan_stim_idx = find(stimpattern(ii,:)>0);
            pulsseq = (1:pulselen);
            stimpattern(ii,chan_stim_idx) = repmat(pulsseq,1,length(chan_stim_idx)/pulselen); % count samples on pulse-pattern
            if ~isempty(chan_stim_idx)
                idx_overlap = chan_stim_idx(ismember(chan_stim_idx,blacklist_idx));
                
                %get complete pulses for stimpattern(idx_overlap) (e.g. if there is only
                % 2 3 in there, get also the 1.
                idx_begin_seq = idx_overlap; %Get starting points of sequences
                idx_overlap_complete = [];
                for dd = 1:length(idx_begin_seq)
                    if pulselen == 4
                        switch stimpattern(ii,idx_begin_seq(dd))
                            case 1
                                idx_overlap_complete(dd,:) = idx_begin_seq(dd):idx_begin_seq(dd)+pulselen-1;
                            case 2
                                idx_overlap_complete(dd,:) = idx_begin_seq(dd)-1:idx_begin_seq(dd)+pulselen-2;
                            case 3
                                idx_overlap_complete(dd,:) = idx_begin_seq(dd)-2:idx_begin_seq(dd)+pulselen-3;
                            case 4
                                idx_overlap_complete(dd,:) = idx_begin_seq(dd)-3:idx_begin_seq(dd);
                        end
                    elseif pulselen == 3
                        switch stimpattern(ii,idx_begin_seq(dd))
                            case 1
                                idx_overlap_complete(dd,:) = idx_begin_seq(dd):idx_begin_seq(dd)+pulselen-1;
                            case 2
                                idx_overlap_complete(dd,:) = idx_begin_seq(dd)-1:idx_begin_seq(dd)+pulselen-2;
                            case 3
                                idx_overlap_complete(dd,:) = idx_begin_seq(dd)-2:idx_begin_seq(dd);
                        end
                    end
                end
                if ~isempty(idx_overlap_complete)
                    idx_overlap = idx_overlap_complete(:);
                end
                stimpattern(ii,idx_overlap) = 0; % set overlapping idx to zero
            end
            blacklist_idx = [blacklist_idx chan_stim_idx]; %Add all newly found pulses to blacklist
            blacklist_idx = sort(unique(blacklist_idx));
        end
        electrodogramm_out = zeros(size(electrodogramm));
        stimpattern(stimpattern>1) = 1; %Remove 2,3s etc.
        electrodogramm_out(logical(stimpattern)) = electrodogramm(logical(stimpattern)); %New electrodogramm
    end

    function [delayed_signal, signal,output_signal] = test_gammatonefilterbankdelay(signal,CISIM_parameter,auralisation_parameter)
        % Create analyse -Filterbank
        analyzer_stim = Gfb_Analyzer_new(CISIM_parameter.voc_sampling_frequency_hz,...
                                            CISIM_parameter.gamma_order_stimulation,...
                                            CISIM_parameter.center_frequencies_hz_stimulation,...
                                            CISIM_parameter.bandwidth_factor);
        % First analysis-filterbank
        %1st analyse:
        [signal, analyzer_stim] = Gfb_Analyzer_process(analyzer_stim,signal); 
        frequency_response = fft(real(signal)');                     
        frequency = [0:length(signal)-1] * auralisation_parameter.voc_sampling_frequency_hz / length(signal);                        
        Gfb_plot(80, [0, auralisation_parameter.voc_sampling_frequency_hz/2, -40, 0], ...
        'frequency response of the analysis filterbank', ...
        'frequency / Hz', 'filter response / dB', ...
        frequency, 20 * log10(abs(frequency_response)));
        hold on; plot((1:length(frequency)),-3.*ones(size(frequency)),'k-');
        disp(' ');
        disp('Figure 80 shows the frequency response of the analysis filterbank.');
        % Create analyse -Filterbank for auralisation
        analyzer_aural = Gfb_Analyzer_new(auralisation_parameter.voc_sampling_frequency_hz,...
                                        auralisation_parameter.gamma_order_auralisation,...
                                        auralisation_parameter.center_frequencies_hz_auralisation,...
                                        auralisation_parameter.bandwidth_factor);
        % Next Analyse-Filterbank (modulates puls-pattern to desired frequency region on basiliar 
        % membrane)
        signal = Gfb_Analyzer_process(analyzer_aural, signal);

        % Signal-Synthesis:
        synthesizer_aural = Gfb_Synthesizer_new (analyzer_aural, 1/100); %1/100 = Delay
        delayed_signal = Gfb_Delay_process(synthesizer_aural.delay, signal);
        % Process signal completely
        output_signal = Gfb_Synthesizer_process(synthesizer_aural, signal);
    end

    function [delay_block] = calculate_delay(fs,pps, n, m)
        factor = m / n; %Adjustment factor for a fixed total stimulation rate
        correct_delay = fs/(pps * factor);
        if correct_delay < 1
            error('Too high pps! This value is not supported!');
        end
        delay_block = ceil(fs/(pps * factor));
        while round(mean(delay_block)*1e5) ~= round(correct_delay*1e5)
            difference = correct_delay - mean(delay_block);
            if difference > 0
                delay_block = [delay_block ceil(fs/(pps * factor))];
            else
                delay_block = [delay_block floor(fs/(pps * factor))];
            end
        end
        if length(delay_block) > 1
            warning('pps * M/N is not a divisor of fs. Make sure that the following equation with x = [1,2,3,4,...,100 is satisfied: x*pps*(M/N) = fs');
        end
    end

    function [pps_hat, ipipps_hat, ipipps_median] = estimate_ppsinterval(electrodogramm, fs)
       % function to estimate the instantaneous pps rate of a given electrodogramm.
       % This is needed to correctly measure T/M-Levels of a patient.
       ipipps_hat = nan(size(electrodogramm));
       ipipps_median = nan(size(electrodogramm,1),1);
       pps_hat = nan(size(electrodogramm,1),1);
       for ii = 1:size(electrodogramm,1)
           tmp_elec = electrodogramm(ii,:);
           % Get all idx of stimulation
           [samplenr] = find(tmp_elec > 0);
           diff_samples = diff(samplenr);
           % remove additional pulse samples
           idx = find(diff_samples == 1); %E.g. adjactend pulsesamples
           tmp_elec(samplenr(idx+1)) = 0; %set additional pulse samples to zero (e.g. no pulse)
           [samplenr] = find(tmp_elec > 0);
           diff_samples = diff(samplenr);
           diff_time = diff_samples./fs; %convert to time values (seconds) 
           % Remove IPIs which are larger than 0.2s Hz, so that our pps
           % estimate is not influenced by pauses in the signal
           silence = find(diff_time>0.2);
%            length_shorted_signal = length(
           pps_hat(ii) = length(samplenr)/length(tmp_elec)*fs;
           ipipps_hat(ii,1:length(diff_time)) = 1./diff_time;
           ipipps_median(ii) = median(ipipps_hat(ii,~isnan(ipipps_hat(ii,:)))); % Calculate median over current electrode pps
       end
    end
%EOF
