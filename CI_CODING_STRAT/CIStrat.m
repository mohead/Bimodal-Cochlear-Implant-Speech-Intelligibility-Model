function [Output, vDistance, StimOrder,vTimeAll, mEnv] = CIStrat (vSignal,par,fs)
% Function to perform standard coding strategies for MED-EL cochlear implants
% and optional Formant locking postprocessing
%
% Input:
% --------
%         vSignal   = Acoustic input sound signal
%         par       = Parameters set in an external parameter file
%         fs        = Sampling frequency of the input wav file
% Output:
% --------
%         Output    = Matrix 12 electrodes x pulses over time [CU]
%         vDistance = Distance between pulse offset and onset of next
%                     pulse [µs]
%         StimOrder = Vector with electrode numbers that are stimulated
%                     over time
%
% Anja Eichenauer, 05.02.2016
% Universität Oldenburg

%% Check Parameters
if par.CIFs~=fs
    vSignal = resample(vSignal,par.CIFs,fs); % resample if fs not equal to CiFs (16000)
end
if size(vSignal,2) > 1
    vSignal=vSignal(:,1); % only monaural
end
if par.m-length(par.nDeactivated)<par.n % sind z.B. 7 Elektroden aktiv, können max 7 aktiviert werden
    par.n=par.m-length(par.nDeactivated);
    error('n larger than available electrodes')
end

if length(par.T) < par.m || length(par.C) < par.m % please set Threshold of inactive Electrode to 0!
    error('Not enough threshold values')
end

if par.T > par.C
    error('Threshold cannot be larger than Comf')
end

if any(par.T < 0) || any(par.C > 1200) % Check if range (0-1200) is eingehalten
    error('T or C Level out of possible range')
end


%% Processing Chain
% Preallocation
FrameIdx=1; % First sample of current frame for block processing
DiffLength=0; %  Length that signal is shorter than frame length
kk=0; % counter
vTimeFrame=0; % center time points of each frame
StimOrder=[];
vDistance=[];
Output=[];
mPulseCycles=[];
vTimeAll=0;
RoundingError=0; % needed for rounding error in Frame shift

%% Processing
FrameLength=0.008; % 0.008 entspricht bei fs=16000 128 Samples
Window=hann(FrameLength*par.CIFs,'periodic');
LastFrame=floor(FrameLength*par.CIFs-min(1./par.pps)*par.CIFs);
while DiffLength <LastFrame % runs until input signal + 127 zeros is over
    kk=kk+1;
    [Input, ~, FrameIdx, DiffLength] = myWindowing (vSignal, par.CIFs, FrameLength,FrameIdx);
    vPreEmp = PreEmphasis(Input); % vPreEmp = emphasized signal frame
    vPreEmp=Window.*vPreEmp; % multiplication with Hann window
    [vSpec, vFreqSpec] = CISpectrum(vPreEmp,par);  % FFT
    [bins, CenterFreq] = CIFilterbank(par);  % Generate Filterbank
    vSigCh = ApplyFilterbankEnv(bins, vSpec); % Apply Filterbank on Signal and calc Env of Channels
    nOfm = SelectCh(vSigCh,par); % Applies n of m procedure
    % Generate stimulation sequences according to parameters
    [PulseCycle, vTime, nSlots, Distance, StimulationOrder]=GeneratePulseCycles(par,nOfm,vTimeAll);
    if par.pps(1) ~= 3000 % 3000 is signiture rate for FL processing
    % Apply loudness-fct. and relate it to T and C-Level
    ElectricSequence =acoustic2electric(PulseCycle, par, nSlots);
    Output(:,size(Output,2)+1:size(Output,2)+nSlots)=ElectricSequence; % after compression
    end
    %% concatenate sequence vectors
%     mPulseCycles(:,size(mPulseCycles,2)+1:size(mPulseCycles,2)+nSlots)=PulseCycle; % before compression
%     vDistance(size(vDistance,2)+1:size(vDistance,2)+nSlots)=Distance'; %
%     StimOrder(1,size(StimOrder,2)+1:size(StimOrder,2)+nSlots)=StimulationOrder;
%     vTimeAll(1,size(vTimeAll,2)+1:size(vTimeAll,2)+(nSlots))=vTime;
    
    %% Frame shift calculation, always reciprocal of stimulation rate
    FrameShift=min(1./par.pps)*par.CIFs; % Periodendauer bestimmt frame shift
    RoundingError=RoundingError+(round(FrameShift)-FrameShift); % accumulation of rounding error
    if abs(RoundingError)>= 1 % because of Rundungsfehler
        FrameShift=(round(FrameShift)-sign(RoundingError));
        RoundingError=RoundingError-sign(RoundingError);
    else
        FrameShift=round(FrameShift);
    end
    if FrameShift > FrameLength*par.CIFs
        error('Frameshift cannot be larger than Framelength')
    end
    FrameIdx=FrameIdx+FrameShift; % Next start frame
    
     %only for plotting purpose
        mSpec(:,kk)= vSpec;
        mEnv(:,kk) = nOfm; %Keeps n-of-m processed envelope, too
end

%% Additional FL Processing
if par.pps(1) == 3000 % 3000 pps == Preprocessed CIS stimulation for FL postprocessing
    FL_ChannelAllocation=2;  % 0 = FL-F0, 1=FL-F1, 2=FL-F0F1F2, 3 = FL-F1-E
    [mFormantData, FormantDataInPl] = ExtractFormantData(vSignal, par,vTimeAll); % Please check Parameters for your sound file! You might have to change window length or nLPC
%     figure;
%     imagesc(vTimeAll, vFreqSpec, 10*log10(max(abs(mSpec).^2,10^(-5))));
%     hold on
%     plot(vTimeAll,FormantDataInPl,'linewidth',3)
%     axis xy

    RateTemplate = GenerateRates(par,mFormantData,FrameShift); 
    [~, NewOutput, NewTime, StimOrder, vDistance,~] = AllocateFormantMod2Ch(CenterFreq, FL_ChannelAllocation, mPulseCycles, RateTemplate, mFormantData,par, nSlots);
 
    vTimeAll=NewTime;
    Output=NewOutput;
end
