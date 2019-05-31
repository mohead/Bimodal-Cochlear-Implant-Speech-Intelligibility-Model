% copy of DemoBSIM, changed for different HRTF, different conditions rigth
% and left
function [SRT] = BimodalSIM(subject_group,varargin)
addpath([fileparts(which(mfilename)) filesep 'Functions']);
addpath([fileparts(which(mfilename)) filesep '..' filesep 'hrir']);
addpath([fileparts(which(mfilename)) filesep '..' filesep 'Vocoder']);
addpath([fileparts(which(mfilename)) filesep '..' filesep 'Vocoder' filesep 'Gammatonefilterbank']);
addpath([fileparts(which(mfilename)) filesep '..' filesep 'openMHA' filesep 'automatic_scripts']);
addpath([fileparts(which(mfilename)) filesep '..' filesep 'openMHA' filesep 'automatic_scripts' filesep 'simu_audibility_hearing_loss']);
addpath([fileparts(which(mfilename)) filesep '..' filesep 'openMHA' filesep 'automatic_scripts' filesep 'simu_distortion_hearing_loss']);


%% RANDOOM NUMBER GENERATOR SEED!!!!
% rng(0) %!!!!!!!!!!!!!!!!!!!!!!!!!!
%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%parse input arguments:
p = inputParser;
p.addRequired('subject_group',@ischar);
p.addParameter('short_time_BSIM_flag',false,@islogical);  %flag to use short-time stBSIM2010 (true) or do Batch Processing BSIM2010 (false)
p.addParameter('noise_azim',[-90,0,90],@isnumeric); % Angle of the noise
p.addParameter('room','anechoic',@ischar); % Room-type, currently only anechoic is supported
p.addParameter('error_flag', true, @islogical);  % flag to use processing errors in EC mechanism (true) or not (false)
p.addParameter('Use_Shadow_filtering', true, @islogical); %flag to use shadow filtering, in MHA preprocessing
p.addParameter('Use_HL_Simulations',false,@islogical); % flag to use hearing loss simulation (true) or not (false)
p.addParameter('plotSRT', false, @islogical); % Plot SRT results
p.addParameter('Display','fig',@(x) (ismember(x,{'text','notext','fig'}))); %Visualization
p.addParameter('plotFigures','none',@(x) (ismember(x,{'none','all','SNR','EC'}))); % 'all': all plots, 'SNR': SNR plots alone, 'EC': EC step plots alone, or 'none' 
p.addParameter('plotFreqBand',8,@(x) (isnumeric(x) && min(x) >= 0 && max(x) <= 30)); % select a band index to plot between 1 and 30, or 0 for all
p.addParameter('plotAngle',-90,@isnumeric);
p.addParameter('Nchanns', 2, @isnumeric);
p.addParameter('Bimodal_SII_Switch_value', 21, @isnumeric);

p.parse(subject_group, varargin{:})
par = p.Results;
   
%% model properties {short-time vs. batch and EC processing error (yes/no)...}

group_para          	= BSIM_subject_group_flags(subject_group,par.Use_HL_Simulations);

EC_flag                 = group_para.EC_flag;   % enable (1) or disable (0) EC processing

better_ear_flag         = group_para.better_ear_flag; % [1=better Ear, 2= Left Ear, 3= Right Ear] (only works for EC_flag = 0)

SII_vals                = group_para.sii_vals;  % SII value. Needs to be adjusted to speech material and SRT of interest

Nchanns                 = 2;  % Numberof Channels: 4channels are needed for ADM algorithm


mfCI_Audiogram          = group_para.Audiogram; % PTAudiogram(group_para.Audiogram);

if isempty(find(par.noise_azim==par.plotAngle,1))
    error(['Please enter a value for plotAngle that belongs to the following set: ' num2str(par.noise_azim)])
end 
%%
% allocate zeros for vector containing SRT values later on.
SRT = zeros(length(par.noise_azim),length(SII_vals));

%% Calibration of signals 
    % calibration of signals to same level for noise and speech (mean between ears)
    ref = 10^(-120/20);                            % corresponds to 120 dB SPL at 0 dB FS
    lev2be = 65;                        % level at which speech is presented
    
    file_path_noise  = ['.' filesep 'olnoise','.wav'];
    file_path_signal = ['.' filesep 'olnoise','.wav'];
     
    [S, fs_S] = audioread(file_path_signal);
    S = S(1:floor(length(S)/2),:);
    S = convolveHRIR(S,fs_S, 0, Nchanns);
    S0 = S;
     
    [N, fs_N] = audioread(file_path_noise);
    N = N(floor(length(N)/2)+1:2*floor(length(N)/2),:);  
    N_ = N;
    N = convolveHRIR(N,fs_N, 0, Nchanns);
      
    lev_S = 20*log10(rms(S)/ref);       % actual rms-level
    lev_S = mean(lev_S);                % frontal speech, should be the same at each ear, reference is MEAN level between the two ears
    Delta_L_speech = lev2be - lev_S;
    
    lev_N = 20*log10(rms(N)/ref);       % actual rms-level
    lev_N = mean(lev_N);                % reference is MEAN level between the two ears
    Delta_L_noise = lev2be - lev_N;  
    % Delta_L values are used for all spatial conditions.
    % Be sure you do NOT calculate Delta_L values for each spatial
    % condition (this would change interaural level differences)
   
for kk=1:length(par.noise_azim)
    
    % Restore Signal to base state with every iteration
    S = S0;
    
    % Load and convolve noise depending on noise_azimuth
    N = N_;     
    N = convolveHRIR(N,fs_N, par.noise_azim(kk),Nchanns);

    %apply calibration
    S = ampSig(S,Delta_L_speech); % amplify each channel accordingly

   
    N = ampSig(N,Delta_L_noise);    
    
    % further processing denpending on subject group
    [S, N] = subject_group_processing(subject_group,S,fs_S,N,fs_N,par.noise_azim(kk),par.Use_Shadow_filtering,par.Use_HL_Simulations);

    % get lengths
    sz_N = size(N);
    sz_S = size(S);
    len_N = max(sz_N);
    len_S = max(sz_S);
    
    % run BSIM: signals should have same length
    len = min([len_N, len_S]);
    N = N(1:len,:);
    S = S(1:len,:);
    
    %apply BSIM level definition
    S = scaleSignalbyRMS(S,120,'dB'); %the level definition of BSIM is 0 dB FS == 0 dB SPL, while the former level definition was == 120 dB SPL
    N = scaleSignalbyRMS(N,120,'dB');
    
    mfSignal = [S, N]; %[SpeechL SpeechR NoiseL NoiseR]
    
    frame_params = [len len/2];
    
    % for short-time stBSIM2010:
    if par.short_time_BSIM_flag
        blocksize = 2.^10; % use the blocksize you prefer (2^10=1024samp=ca.23ms @ 44100 Hz)
        frame_params = [blocksize blocksize/2];
    end
    
    % processing
    [SRT(kk,1:length(SII_vals)),stData,frameData] = runBSIM(mfSignal,'frame',frame_params,'display',par.Display,'errorflag',par.error_flag,'CI_audiogram',mfCI_Audiogram,'ECprocessing',EC_flag,'BetterEar',better_ear_flag,'angle',par.noise_azim(kk),'short_term',par.short_time_BSIM_flag,'SIIparam',SII_vals, 'subject_group', subject_group,'plotfigures',par.plotFigures,'plotfreqband',par.plotFreqBand,'plotangle',par.plotAngle,'Bimodal_SII_Switch_value',par.Bimodal_SII_Switch_value);
end

% plot the results in the way it is done in the article:
if par.plotSRT
    plotSRT_angle (par.noise_azim,SRT,[par.room '  ' subject_group],mfCI_Audiogram,subject_group,SII_vals)
%     savefig(['..' filesep 'Figures' filesep subject_group '_SRT.fig']);
end 
end
% Plotting
function plotSRT_angle (angles,SRT,TitleText,mfAudiogram, subject_group,SII_vals)

figure('Name', ['''' subject_group '''']);
subplot(1,2,1);
plot(angles,SRT,'-o','Color','k','Linewidth',2,'Markersize',9);
set(gca,'XTick',[-180 -90 0 90 180]);
set(gca,'YTick',-27:3:144);
ylim([-24 9]);
xlim([-185 185])
grid on;
title(TitleText);
xlabel('\bf{Noise azimuth / degree}');
ylabel('\bf{SRT / dB SNR}');
text(-170, 8,['SII ' subject_group ' = '],'HorizontalAlignment','left');

text(-170, 5,num2str(SII_vals),'HorizontalAlignment','left');
subplot (1,2,2);
plotaudi(mfAudiogram(:,1), mfAudiogram(:,2:3), 'lr', TitleText);

end