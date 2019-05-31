function [output_sig_hg] = hearing_loss_sim(signal,fs,HGType,phi,Use_Shadow_filtering,Use_HL_Simulations,Signal_type)
% Simulation of hearing loss




freq_hl_params = freq_distortion_params;
freq_hl_params.numIterations = 2;
freq_hl_params.bandwidthFactor = 2; % Bandwidth of HI filterbank is just 2 instead of 4 (default)
if Use_Shadow_filtering
    
if phi==0
    phi_string = '000_';
elseif phi==90
    phi_string = 'p90_';
elseif phi==-90
    phi_string = 'm90_';
end

if strcmp( HGType,'MHL')
    HL_string = 'HL2';
elseif strcmp( HGType,'SHL')
    HL_string = 'HL1';
end 

output_sig_hg      	= audioread([Signal_type phi_string HL_string '.wav']);

else
% Obtain current directory
    mfile_full_path = mfilename('fullpath');
    OPENMHA_dir     = fileparts(fileparts(mfile_full_path));
    switch HGType
        case 'SHL'
            mha_config_file    = [OPENMHA_dir filesep 'automatic_scripts' filesep 'HG1_fitting_CAMFIT_compromise_openMHA.cfg'];
        case 'MHL'
            mha_config_file    = [OPENMHA_dir filesep 'automatic_scripts' filesep 'HG2_fitting_CAMFIT_compromise_openMHA.cfg'];
        case 'NH'
            mha_config_file    = [OPENMHA_dir filesep 'automatic_scripts' filesep 'NH_nogain.cfg'];
        case 'test'
            mha_config_file    = [OPENMHA_dir filesep 'automatic_scripts' filesep 'HG2_fitting_CAMFIT_compromise_openMHA.cfg'];
        otherwise
            error('HGType must be one of the following ''HG1'',''HG2''');
    end
    addpath([OPENMHA_dir filesep 'matlab']);
    audioFolder = [OPENMHA_dir filesep 'mha' filesep];
    infile = strcat(audioFolder,'MHA-tmp-in.wav');
    outfile = strcat(audioFolder,'MHA-tmp-out.wav');
    % scale input signal to MHA_calib_input (120 dB SPL)
    %signal = scaleSignalbyRMS(signal,-mha_calib_correction_factor,'dB');
    %write signal to disk:
    audiowrite(infile, signal, fs, 'BitsPerSample', 32); %Input to MHA needs to have 4 channels, even when applying no processing
    % run MHA processing
    % worders_ID = parpool(4);
    call_mha_generic('../openMHA/mha/configurations/dynamiccompression.cfg',infile,outfile,'mha_dir','/home/ayham/Code/BSc-Thesis_Anja_Kalin/openMHA/bin/', 'plugin_name','mha.overlapadd.mhachain.dc','plugin_config_file',mha_config_file);
    % delete worders_ID
    % read in resulting signal
    [mha_sig,fs] = audioread(outfile);
    % Return Signal to BSIM level
    %mha_sig = scaleSignalbyRMS(mha_sig,mha_calib_correction_factor,'dB'); % (Assumed: MHA_calib_output is the same as input (120 dB SPL)
    output_sig_hg = mha_sig;

end
if Use_HL_Simulations %only perform impairment simulation if requested   
    %% Simulation processing
    switch HGType
        case 'SHL'
            hearing_thresholds_oet = [45 70 100 100 100 100]; % dB HL for 250, 500, 1000, 2000, 4000, 8000 Hz
        case 'MHL'
            hearing_thresholds_oet = [0 0 30 50 60 70]; % dB HL for 250, 500, 1000, 2000, 4000, 8000 Hz
        case 'NH'
        case 'test'
            hearing_thresholds_oet = [0 0 30 50 60 70]; % dB HL for 250, 500, 1000, 2000, 4000, 8000 Hz
        otherwise
            error('HGType must be one of the following ''HG1'',''HG2''');
    end
    % if def.simulate_hearing_loss
    output_sig_dist = f_process(output_sig_hg(:,1),freq_hl_params); %Frequency distortion component by Robert Bennett-Algorithm
    output_sig_dist = output_sig_dist.synthesis(:); %This is the frequency distorted signal
    % debug('Finished frequency distortion simulation',param,@analyse_in_out,signal,output_sig_dist(1:length(signal)),fs, 'Original (1) vs. Frequency smearing hearing loss (2)')

    % add the db value for the oetting algo (Goal of 65 dB SPL FF, e.g.
    % rms2(signal,1) = 65 dB for a reference signal (olnoise, from the front)
    output_sig_dist = scaleSignalbyRMS(output_sig_dist,120+3,'dB'); % +2 because to account for bennett distortion hearing simulation level loss.
    output_sig_hg   = simulate_hearing_loss(output_sig_dist,fs,hearing_thresholds_oet,freq_hl_params.debug); %Hearing loss audibility component
    % subtract the diff dB value from the signal again
    output_sig_hg   = scaleSignalbyRMS(output_sig_hg,-120,'dB');
    % output_sig_dist = scaleSignalbyRMS(output_sig_dist,-oetting_diff_dB,'dB');

    % debug('Finished acoustic audibility simulation',param,@analyse_in_out,output_sig_dist,output_sig_hg(1:length(signal)),fs, 'Frequency smearing hearing loss (1) vs. Frequency smearing hearing loss + Audibility hearing loss (2)')
end
 
% % %% Simulation processing
% % % To be added soon.
end
