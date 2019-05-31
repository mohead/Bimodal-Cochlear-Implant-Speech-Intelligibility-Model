function [Output, OutputFL, NewTime, StimOrder, vDistance,AllFormantCh2] = AllocateFormantMod2Ch(CenterFrqs, FL_ChannelAllocation,EnvSig,FormantTimingInter, mFormantData, par, nSlots)
% Function to allocate channels according to F0 and Formant frequencies
%
% Input:
% --------
%         CenterFrqs            = Center frequs of filter bank
%         FL_ChannelAllocation  = 0 = FL-F0, 1=FL-F1, 2=FL-F0F1F2, 3 = FL-F1-E
%         EnvSig                = Pure Electrodogram before compression
%         FormantTimingInter    = Rate template   
%         mFormantData          = [Time; VoicedUnvoiced; F0, F1; F2] interpolated
%         par                   = Parameters set in an external parameter file
%         nSlots                = no of slots per cycle
% Output:
% --------
%         Output                = Matrix 12 electrodes x pulses over time [CU]
%         OutputFL              = Electrodogram after FL processing
%         NewTime               = Time vector after cutting procedure 
%         StimOrder             = Vector with electrode numbers that are stimulated over time
%         vDistance             = Distance between pulse offset and onset of next pulse [µs]
%         AllFormantCh2         = First Formant related channel (transition channels)
%
% Anja Eichenauer, 05.02.2016
% Universität Oldenburg

Time=mFormantData(1,:);
Count=0;
NewTime=0;
DistanceDifference=0;
Output=zeros(size(EnvSig));
OutputFL=zeros(size(EnvSig));
VoicedUnvoiced=mFormantData(2,:);
F0InPl=mFormantData(3,:);
F1InPl=mFormantData(4,:);
F2InPl=mFormantData(5,:);
AllFormantCh=zeros(size(F1InPl,2),3);
% Allocate Channel to Formant frequency
for FChannels=2:size(EnvSig,2)
    if VoicedUnvoiced(FChannels)==0 % Unvoiced
        F0Ch=1; % Channels don´t make a difference here,
        F1Ch=2; % ... ...because if not voiced, all channels same rate.
        F2Ch=3;
        AllFormantCh(FChannels,:)=[F0Ch F1Ch F2Ch];
        FormantFrequencies=900; % only standard CIS
    else
        [~, F0Ch]=min(abs(CenterFrqs-F0InPl(FChannels))); % find channel closest to desired frequency
        [~, F1Ch]=min(abs(CenterFrqs-F1InPl(FChannels)));
        [~, F2Ch]=min(abs(CenterFrqs-F2InPl(FChannels)));
        AllFormantCh(FChannels,:)=[F0Ch F1Ch F2Ch];
        FormantFrequencies=[F0InPl(FChannels) F1InPl(FChannels) F2InPl(FChannels)];
        if F0Ch >= F1Ch
            F1Ch=F0Ch+1;
        end
        if F1Ch>= F2Ch
            F2Ch=F1Ch+1;
        end
    end
    
    % matrix of different rate-templates is generated according to
    % electrode-frequency allocation:
    switch FL_ChannelAllocation %0 = FL-F0, 1=FL-F1, 2=FL-F0F1F2, 3 = FL-F1-E
        % How many Formants are Coded
        case 0
            nF0=3; % no of channels for F0-coding
            FixedChannelAllocation=[FormantTimingInter(1,FChannels)*ones(nF0,1);  FormantTimingInter(3,FChannels)*ones(12-nF0,1)];
            FormantChannels=nF0;
        case 1
            nF0=0;
            nF1=12; % All channels F1 pps
            FixedChannelAllocation=[FormantTimingInter(1,FChannels)*ones(nF0,1);...
                FormantTimingInter(2,FChannels)*ones(nF1,1)];
            FormantChannels=[nF0 nF1];
        case 2
            nF0=F1Ch-1; % F0 always starts at Ch1
            nF1=F2Ch-F1Ch; %number of channels for each rate
            nF2=13-F2Ch;
            FixedChannelAllocation=[FormantTimingInter(1,FChannels)*ones(nF0,1);...
                FormantTimingInter(2,FChannels)*ones(nF1,1); FormantTimingInter(3,FChannels)*ones(nF2,1)];
            FormantChannels=[nF0 nF1 nF2];
        case 3
            nFX=7; % 7 channels 900 pps fix
            nF1=5; % 5 channels F1-E pps
            FixedChannelAllocation=[FormantTimingInter(2,FChannels)*ones(nF1,1); FormantTimingInter(3,FChannels)*ones(nFX,1)];
            FormantChannels=[0 nF1 nFX];
    end
    % Adaptive T and C value allocation
    ElectricSequence =acoustic2electricAdaptive(EnvSig(:,FChannels), par, nSlots, FormantChannels, FormantFrequencies); % New: Thresholds are adapted here according to FormantFreq
    Output(:,FChannels)=ElectricSequence.*FixedChannelAllocation; % Amplituden das Muster aufprägen

    %% Pulse shifting
    % "Cutting" Procedure, Signal is sampled with 3000 pps, then only
    % wanted pulses are cut out of this Electrodogram
    % Necessary to be able to use 40us pulse width
    MinDist=(2*par.PulseWidthWrapper+(8*10^-6));
    [row ~]=find(EnvSig(:,FChannels)); % which channel is active?
    if any(row)==0
        row=1;
    end
    if FixedChannelAllocation(row)==1 % Is there a pulse at this time?
        Count=Count+1;
        NewTime(Count+1)=Time(FChannels); % Only save time values when pulses occur
        AllFormantCh2(Count,:)=AllFormantCh(FChannels,:);
        StimOrder(Count)=row; % update stim order
        OutputFL(row,Count)=Output(row,FChannels); % update electrodogram
        vDistance(Count)=diff(NewTime(end-1:end)); % update distance
        if vDistance(end)<=MinDist % minimal distance necessary because of pulsewidth (40 us*2+ipg + puffer)
            DistanceDifference(Count)=(abs(diff([MinDist vDistance(end)]))); % sum up adjustments
            NewTime(Count+1)=NewTime(Count+1)+DistanceDifference(Count); % shift time vector
            vDistance(end)=MinDist; % new distance
        end
        
    end
end

ElectrodogramEnd=size(StimOrder,2); % Cut end becuase it is shorter than the original
OutputFL=OutputFL(:,1:ElectrodogramEnd);
NewTime=NewTime(2:ElectrodogramEnd(end)+1);