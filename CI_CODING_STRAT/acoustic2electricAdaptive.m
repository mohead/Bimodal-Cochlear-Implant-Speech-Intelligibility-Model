function Levels =acoustic2electricAdaptive(PulseCycle,par,~,FormantChannels, FormantFrequencies)
% Function to allocate channels according to F0 and Formant frequencies
%
% Input:
% --------
%         PulseCycle            = Electrodogram
%         par                   = Parameters set in an external parameter file
%         FormantChannels       = Allocated Electrodes
%         FormantFrequencies    = Formant frequencies

% Output:
% --------
%         Levels                = Current units for all 12 electrodes [CU]

%
% Anja Eichenauer, 05.02.2016
% Universität Oldenburg

%% Threshold allocation
T=zeros(12,1); % only stimulation with 900 pps
C=zeros(12,1);
LoudnessFactor=par.vol;
T250=par.T(1:12)*LoudnessFactor; % Thresholds 250
T400=par.T(13:24)*LoudnessFactor; % Thresholds 500
T900=par.T(25:36)*LoudnessFactor;  % Thresholds 900
C250=par.C(1:12)*LoudnessFactor; % Thresholds 250
C400=par.C(13:24)*LoudnessFactor; % Thresholds 500
C900=par.C(25:36)*LoudnessFactor;  % Thresholds 900
if size(FormantFrequencies,2) == 1 &&  FormantFrequencies == 900 % in unvoiced case
    T(1:12)=T900; % only stimulation with 900 pps
    C(1:12)=C900;
else
    T(1:FormantChannels(1))=T250(1:FormantChannels(1)); % always 250-Schwelle for F0
    C(1:FormantChannels(1))=C250(1:FormantChannels(1)); % always 250-Schwelle for F0
    if size(FormantChannels,2)==1 % wenn nur F0
        T(FormantChannels(1)+1:12)=T900(FormantChannels(1)+1:12); % Rest mit 900 pps-Schwelle
        C(FormantChannels(1)+1:12)=C900(FormantChannels(1)+1:12); %
    else % F1
        if FormantFrequencies(2) <= 500 % Threshold E I O U
            T(FormantChannels(1)+1:FormantChannels(1)+FormantChannels(2))=T400(FormantChannels(1)+1:FormantChannels(1)+FormantChannels(2)); % Threshold for low frequency F1
            C(FormantChannels(1)+1:FormantChannels(1)+FormantChannels(2))=C400(FormantChannels(1)+1:FormantChannels(1)+FormantChannels(2)); %
            if size(FormantChannels,2)==3
                T(FormantChannels(1)+FormantChannels(2)+1:12)=T900(FormantChannels(1)+FormantChannels(2)+1:12); % F2 Channels are stimulated with 900 pps
                C(FormantChannels(1)+FormantChannels(2)+1:12)=C900(FormantChannels(1)+FormantChannels(2)+1:12); %               
            end
        else % Threshold A, if Formant Freq > 500
            T(FormantChannels(1)+1:12)=T900(FormantChannels(1)+1:12); % Schwelle für A und Rest gleich
            C(FormantChannels(1)+1:12)=C900(FormantChannels(1)+1:12); %
        end
    end
    
end

%%
NoSigElectrodes=PulseCycle<par.B; % Electrodes that do not contain any audible Signal
PulseCycle(PulseCycle< par.B)=par.B; % Adjust Amplitudes according to minimum....
PulseCycle(PulseCycle> par.M)=par.M; %... and maximum Amp (Clipping).
CompressedEnv=log(1+par.CompFac*((PulseCycle-par.B)/(par.M-par.B)))/log(1+par.CompFac); % loudnes function
TCLevel=(C-T); % part 1 formular 6
Levels=T+(TCLevel.*(CompressedEnv)); % part 2 formular 6
Levels(NoSigElectrodes)=0; % this is needed because else low signal increments would be set to the threshold level
