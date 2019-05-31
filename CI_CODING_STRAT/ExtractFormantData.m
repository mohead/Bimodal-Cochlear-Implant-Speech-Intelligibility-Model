function [mFormantData, FormantDataInPl] = ExtractFormantData(vSignal, par, vTimeAll)
% Function to extract F0 and Formants (F1 and F2)
%
% Input:  
% --------
%         vSignal   = Acoustic input sound signal
%         par       = Parameters set in an external parameter file
%         vTimeAll  = time vector of CIS stimulation
% Output: 
% --------
%         mFormantData    = [Time; VoicedUnvoiced; F0, F1; F2] interpolated
%         FormantDataInPl = interpolated [F0; F1; F2];
%                    
%
% Anja Eichenauer, 20.02.2016
% Universität Oldenburg

%% Preallocation
alpha=0.3; % recursive smoothing factor
FrameL=0.032; % [s] 
FrameS=FrameL*0.4; % 40 % of frame length

FrameShift=floor(FrameS*par.CIFs);   % FrameShift s to samples
NumOfFrames=ceil(length(vSignal)/FrameShift); % Calculate Nr. of Frames
FrameLength=FrameL*par.CIFs; % transform FrameLength from s to samples

%% Initialise
f0=zeros(NumOfFrames,1);
Formants=zeros(NumOfFrames,2);
sFormants=zeros(1,2); % smoothed
sf0=0; % smoothed

FrameIdx=1; % Start Sample of Frame
mFrames=zeros(FrameLength,NumOfFrames);
vTimeFormant=zeros(1,NumOfFrames);
MeanFormants=[];
FormerFormants=[];
NVoiced=0;
%% Windowing, Voiced/unvoiced, F0, Formant Estimation and smoothing
for jj=1:NumOfFrames
    [mFrames(:,jj), vTimeFormant(jj), FrameIdx, DiffLength] = Windowing (vSignal, par.CIFs, FrameLength,FrameIdx);
    FrameIdx=FrameIdx+FrameShift;
    [VoicedUnvoiced(jj), NZero(jj)] = VoicedUnvoicedSegment(mFrames(:,jj)); % voiced/ unvoiced estimation via zero crossings
    if VoicedUnvoiced(jj)== 1 % if voiced
        f0(jj)= fZeroEstimator(par.CIFs, mFrames(:,jj)); % estimate f0 via autocorrelation
        [Formants(jj,:)]= FormantEstimator(par.CIFs, mFrames(:,jj)); % estimate formants via LPC
        % smoothing
        if jj > 1
            sFormants(jj,:)=(alpha).*sFormants(jj-1,:) + (1-alpha).*Formants(jj,:); % apply recursive smoothing, for very short stimuli better no smoothing(alpha=0)
            sf0(jj)=(alpha).*sf0(jj-1) + (1-alpha).*f0(jj);
        else
            sFormants(jj,:)=Formants(jj,:); % First value without smoothing
            sf0(jj)=f0(jj);
        end
    else % if unvoiced
        f0(jj)= 100; % value only set to fill sample, only needed for very first sample
        Formants(jj,:)= ones(1,2)*200; % value only set to fill sample, only needed for very first sample
        if jj > 1
            sFormants(jj,:)=sFormants(jj-1,:); % no update (smoothing), take previous value
            sf0(jj)=sf0(jj-1);
        else
            sFormants(jj,:)=Formants(jj,:); % First value
            sf0(jj)=f0(jj);
        end
        
    end

end

VoicedUnvoiced= interp1(vTimeFormant,VoicedUnvoiced,vTimeAll,'pchip','extrap');
VoicedUnvoiced=round(VoicedUnvoiced); % Eliminate error from interpolation
F0InPl= interp1(vTimeFormant,sf0,vTimeAll,'linear','extrap');
F0InPl(F0InPl<0)=0;
F1InPl= interp1(vTimeFormant,sFormants(:,1),vTimeAll,'linear','extrap');
F1InPl(F1InPl<0)=0;
F2InPl= interp1(vTimeFormant,sFormants(:,2),vTimeAll,'linear','extrap');
F2InPl(F2InPl<0)=0;
FormantDataInPl=[F0InPl; F1InPl; F2InPl];

mFormantData=[vTimeAll;VoicedUnvoiced; F0InPl; F1InPl; F2InPl]; % matrix with retreived Data for further processing