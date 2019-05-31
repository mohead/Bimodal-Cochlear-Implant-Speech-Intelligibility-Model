function Levels =acoustic2electric(PulseCycle,par,nSlots)
%% function to apply loudness function and map envelope to individual dynamic
% INPUT
% ---------
%   PulseCycle   :    Vector with Amplitudes of this cycle [CU]
%   par          :    Parameters from file CIParameters
%   nSlots       :    Equals n
%   
%
% OUTPUT
% ---------
%  Levels             :   Vector with individual Amplitudes [CU]
if ~isempty(PulseCycle)
NoSigElectrodes=PulseCycle<par.B; % Electrodes that are below base level
PulseCycle(PulseCycle< par.B)=par.B; % Adjust Amplitudes according to minimum....
PulseCycle(PulseCycle> par.M)=par.M; %... and maximum Amp (Clipping).
CompressedEnv=log(1+par.CompFac*((PulseCycle-par.B)/(par.M-par.B)))/log(1+par.CompFac); % loudnes growth function
CompressedEnv=CompressedEnv.*par.vol;
TCLevel=(par.C-par.T); % part 1 formular 6 (Nogueira), equals dynamic range
TCLevel=repmat(TCLevel,1,nSlots); % adjust matrix size to be able to multiply in next step
TLevel=repmat(par.T,1,nSlots); % adjust matrix size to be able to multiply in next step
Levels=TLevel+(TCLevel.*(CompressedEnv)); % part 2 formular 6
Levels(NoSigElectrodes)=0; % this is needed because else signal increments below base would be set to the threshold level
else
    Levels = double.empty(1,0); %Return nothing
end