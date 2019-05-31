function [sOutput] = f_process(in, parameters)

% Author: R. Bennett, 2011, modified by B. Williges, 2016

in = in(:); %Make vector
	
delay_samples = round(parameters.filterDelaySec.*parameters.sampleFreqHz);

%in = in./max(abs(in));

in = [in; zeros(sum(delay_samples),1)]';
%in = hilbert(in);
debug('Input signal',parameters,@plotChannels,real(in),parameters.sampleFreqHz,1,'r','Input signal');
NH_in = in;
HI_in = in;

%% create Filterbanks
% Analysefilterbank objekt Normalhörend (NH)
analyzer_NH = Gfb_Analyzer_new_HLS(parameters.sampleFreqHz,         ...
                                     parameters.lowerCutoffFreqHz,     ...
                                     parameters.baseFreqHz, ...
                                     parameters.upperCutoffFreqHz,     ...
                                     parameters.filtersPerERBaud,           ...
                                     parameters.gammaOrderNH12,             ...
                                     parameters.bandwidthFactor_NH);       ...
                                     %normalization_divisor
analyzer_NH.fast = 0; %Use fast c++ implementation
% NH-Analysefilterbank 1 und 2
[analyzer_NH_1, analyzer_NH_2] = ...
    Gfb_Analyzer_split(analyzer_NH, parameters.splitOrder1, parameters.splitOrder2 );

% Code supplied by Volker Hohmann and addapted by RB start

		delay1 = Gfb_Delay_new(analyzer_NH_1,delay_samples(1));
		delay1_2 = Gfb_Delay_new(analyzer_NH,sum(delay_samples));
		delay2 = delay1; % only init
		delay2.delays_samples = delay1_2.delays_samples - delay1.delays_samples;
		if min(delay2.delays_samples) < 0
		    disp('Delays smaller zero, increase total delay!')
		end
		delay2.memory = zeros(length(delay2.delays_samples), ...
		                      max(delay2.delays_samples));
		delay2.phase_factors = delay1_2.phase_factors ./ delay1.phase_factors;
% Code supplied by Volker Hohmann and addapted by RB stop
 
% Analysisfilterbank hearing impaired (HI)
analyzer_HI = Gfb_Analyzer_new_HLS(parameters.sampleFreqHz,         ...
                                     parameters.lowerCutoffFreqHz, ...
                                     parameters.baseFreqHz, 	...
                                     parameters.upperCutoffFreqHz, ...
                                     parameters.filtersPerERBaud,  ...
                                     parameters.gammaOrderHI,      ...
                                     parameters.bandwidthFactor);  ...
                                     % normalization_divisor
analyzer_HI.fast = 0; %User fast c++ implementation
% Delay HI is the same length as the delay for NH
delay_HI = Gfb_Delay_new(analyzer_HI,delay_samples(1));

mixer_NH = Gfb_Mixer_new_HLS(analyzer_NH, delay1_2);

%% Processing of signals------------------------------------------
% Analyze Signal HI hearing Complete
[output_HI, analyzer_HI,delay_HI] = process_Gammatonefilter(analyzer_HI, delay_HI, HI_in);
debug('Processing hearing impaired output after hearing impaired filterbank',parameters,@plotChannels,output_HI,parameters.sampleFreqHz, 1, 'Processing hearing impaired output after hearing impaired delay');

% Calculate Envelope HI
[Env_HI,phase_HI,TFS_HI] = extract_Env_TFS(output_HI);
    debug('Extracted HI filter output: Temporal FineStructure',parameters,@plotChannels,TFS_HI,parameters.sampleFreqHz,0.45,'Extracted HI filter output: Temporal FineStructure (scaled by 0.45)');
    debug('Extracted HI filter output: Envelope',parameters,@plotChannels,Env_HI,parameters.sampleFreqHz,1,'Extracted HI filter output: Envelope');

for iter = 1:parameters.numIterations
	% ------------------------------------------------------------
	% Analyze Signal NH hearing first split stage
    [output_NH, analyzer_NH_1,delay1] = process_Gammatonefilter(analyzer_NH_1, delay1, NH_in);
	debug('Processed normal hearing output after first NH filterbank',parameters,@plotChannels,output_NH,parameters.sampleFreqHz, 1, 'Processing normal hearing output after hearing impaired delay');
    
    % Get RMS after first NH filterbank in order to not influence the
    % overall level by the number of iterations
    if iter == 1
        rms_ref = rms2(real(output_NH),2);
    end
	% ------------------------------------------------------------
	% Calculate Envelope and Finestructure for NH
	% The Envelope is used to calculate the correlation coefficient
    [Env_NH,~,TFS_NH] = extract_Env_TFS(output_NH);
    debug('Extracted NH filter output: Temporal FineStructure',parameters,@plotChannels,TFS_NH,parameters.sampleFreqHz,0.45,'Extracted NH filter output: Temporal FineStructure (scaled by 0.45)');
    debug('Extracted NH filter output: Envelope',parameters,@plotChannels,Env_NH,parameters.sampleFreqHz,1,'Extracted NH filter output: Envelope');
	% ------------------------------------------------------------
	% Construct 'Chimaera' / Process Stage
	Chimaera = TFS_NH.*Env_HI;
	debug('Constructed Chimaera: TFS_NH * Env_HI',parameters,@plotChannels,Chimaera,parameters.sampleFreqHz,1,'Constructed Chimaera: TFS_NH * Env_HI');

	% ------------------------------------------------------------
	% Filter 'Chimaera' with NH-Filter, second split section to suppress unwanted 
	% sidebands produced by modulation with envelope of HI filter bank
    [Chimaera, analyzer_NH_2,delay2] = process_Gammatonefilter(analyzer_NH_2, delay2, Chimaera);
    debug('Processed Chimaera after second NH filterbank',parameters,@plotChannels,Chimaera,parameters.sampleFreqHz,1,'Processed Chimaera with second NH filterbank delay');
	%-------------------------------------------------------------

	% Set rms in each filter to rms of first NH filterbank
    Chimaera = setSignaltoRMS(Chimaera,rms_ref,'linear');
	
	% Mixer Object with Complex output. For real output use ...mixer2
	% Complex output required here
	[output, mixer_NH] = Gfb_Mixer_process_2(mixer_NH, Chimaera);
    debug('Processed Chimaera after NH synthesis filterbank',parameters,@plotChannels,output,parameters.sampleFreqHz,1,'Processed Chimaera with NH synthesis filterbank');
	NH_in = [output(sum(delay_samples):end), zeros(1,sum(delay_samples)-1)];
	
	% Aligne output signal with input signal to realine at input of following iteration
    [NH_in] = align_signals(in, NH_in); NH_in = NH_in.'; %Nonconjugate transpose (.') is important, otherwise we have a sign shift from + to - in the complex part!
	debug('Delayed NH input for next iteration',parameters,@plotChannels,output,parameters.sampleFreqHz,1,'Delayed NH input for next iteration');

% 	sOutput.iteration(iter).synthesis = NH_in; % output before reanalysis

    fprintf('Completed iteration No %d\n', iter)
end %loop

fprintf('\n')

sOutput.in = in;
sOutput.synthesis = NH_in; %This is the frequency distorted signal
sOutput.output_HI = output_HI; % Hearing impaired analysis output
if parameters.debug
sOutput.corr_coeffs = verify_algo(in,NH_in,parameters); %This does the re-analyis and calculates the cross correlation across the output channels.
end
sOutput.center_frequ_hz = analyzer_NH.center_frequencies_hz;

%% verify function for reanalysis of the simulation output and cross correlation calculation afterwards
    function [cross_correlation] = verify_algo(orig_sig,output_algo, parameters)
       % create NH and SH filterbank
       final_analyzer_NH = Gfb_Analyzer_new_HLS(parameters.sampleFreqHz,      ...
                                     parameters.lowerCutoffFreqHz, ...
                                     parameters.baseFreqHz, ...
                                     parameters.upperCutoffFreqHz, ...
                                     parameters.filtersPerERBaud,  ...
                                     parameters.gammaOrderNH,      ...
                                     parameters.bandwidthFactor_NH);
        final_analyzer_SH = Gfb_Analyzer_new_HLS(parameters.sampleFreqHz,         ...
                                     parameters.lowerCutoffFreqHz, ...
                                     parameters.baseFreqHz, 	...
                                     parameters.upperCutoffFreqHz, ...
                                     parameters.filtersPerERBaud,  ...
                                     parameters.gammaOrderHI,      ...
                                     parameters.bandwidthFactor);
        % Delay for re_analysis, to obtain same number delay samples as HI path
        delay_in_samples = round(parameters.filterDelaySec.*parameters.sampleFreqHz);
        final_delay_NH = Gfb_Delay_new(final_analyzer_NH,delay_in_samples(1)); 
        final_delay_SH = Gfb_Delay_new(final_analyzer_SH,delay_in_samples(1)); 
       
       % process input through both filterbanks
       [reanalyzed_output_NH, final_analyzer_NH, final_delay_NH] = process_Gammatonefilter(final_analyzer_NH, final_delay_NH, output_algo);
       [reanalyzed_output_SH, final_analyzer_SH, final_delay_SH] = process_Gammatonefilter(final_analyzer_SH, final_delay_SH, orig_sig);
       debug('Processed signal after NH reanalyis filterbank',parameters,@plotChannels,reanalyzed_output_NH,parameters.sampleFreqHz,1,'Processed signal after NH reanalyis filterbank');
       debug('Processed signal after SH reanalyis filterbank',parameters,@plotChannels,reanalyzed_output_SH,parameters.sampleFreqHz,1,'Processed signal after SH reanalyis filterbank');
       % calculate Envelope and TFS from both signals
       [Env_NH_re,phase_NH_re,TFS_NH_re] = extract_Env_TFS(reanalyzed_output_NH);
       [Env_SH_re,phase_SH_re,TFS_SH_re] = extract_Env_TFS(reanalyzed_output_SH);
       % calculate cross correlation for the output
       [cross_correlation] = calculate_correlation(Env_NH_re,Env_SH_re);
    end

%% Support functions
    function [output_sig, analyzer,delay] = process_Gammatonefilter(analyzer, delay, inputSignal)
        [output_sig, analyzer] = Gfb_Analyzer_process(analyzer, inputSignal);
        analyzer = Gfb_Analyzer_clear_state(analyzer); % clear state for analyzer
        % Output_HI_delay is the signal that will be used for the iteration process
        [output_sig, delay] = Gfb_Delay_process_2(delay, output_sig); %Old code Gfb_Delay_process_2 (uses complex values for multiplication with phase factor instead of real values)
        delay = Gfb_Delay_clear_state(delay); % clear state for delay object
    end

    function [Env,phase,TFS] = extract_Env_TFS(gammatone_output)
            assert(isreal(gammatone_output) == 0,'Input must be complex');
            Env = abs(gammatone_output);
            phase = angle(gammatone_output);
            TFS = cos(phase);% Reconstructed real signal is then Env.*TFS = real(gammatone_output) (with rounding errors ~10^-16)
    end

    function [corr_coeffs] = calculate_correlation(sig_1,sig_2)
	size_of_channels = min(size(sig_1));
    corr_coeffs_temp = zeros(size_of_channels,2,2);
    corr_coeffs = zeros(size_of_channels,1);
	for channels = 1:size_of_channels
        % calculate time-shift between signals and adjust for it.
        [sig_2_shift] = align_signals(sig_1(channels,:), sig_2(channels,:));

        % Calculation of the correlation coefficient for each channel
        % Shifted by xcorr lags
        corr_coeffs_temp(channels,:,:) = corrcoef((sig_1(channels,:))', sig_2_shift');
        corr_coeffs(channels) = corr_coeffs_temp(channels,2,1);
	end 
    end

    function [aligned_signal] = align_signals(shifted_signal, ref_signal)
        shifted_signal = shifted_signal(:); %Make vector
        ref_signal = ref_signal(:); %Make vector
        [xcorrelation, xlags] = xcorr(real(shifted_signal),real(ref_signal));
        [~, xInd] = max(xcorrelation);
        xlags_in_out = xlags(xInd);
        % adjust ref signal accordingly
            if xlags_in_out < 0
                aligned_signal = [ref_signal(abs(xlags_in_out)+1:end); zeros(abs(xlags_in_out),1)];
            elseif xlags_in_out > 0
                aligned_signal = [zeros(abs(xlags_in_out),1); ref_signal(1:end-abs(xlags_in_out))];
            elseif xlags_in_out == 0
                aligned_signal = ref_signal;
            end
        aligned_signal = aligned_signal(:); %make vector
    end
end