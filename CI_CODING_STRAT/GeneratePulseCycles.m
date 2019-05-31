function [ElectrodeScheme, vTimeSlots, nSlots, Distance, StimulationOrder] = GeneratePulseCycles(par,nOfm,vTimeOld)
%% function to generate stimulation sequence for one signal frame
% INPUT
% ---------
%   par     :    Parameters from file CIParameters
%   nOfm    :    Vector with Amplitudes of this cycle
%   vTimeOld:    Time vector of previous cycle [s]
%
% OUTPUT
% ---------
%  ElectrodeScheme    :   Matrix 12 Electrodes x n pulses
%  vTimeSlots         :   Time vector of current cycle [s]
%  nSlots             :   Equals n
%  Distance           :   Vector with time between 2 pulses [s]
%  StimulationOrder   :   Vector with Electrode numbers in stimulation
%                         order

PulsePeriod=1./par.pps(1);
Distance=(PulsePeriod/par.n); % Distance between two pulses
if Distance <= par.PulseWidth*10^-6
     error('Pulse Width too large for this stimulation rate')
end
% Stimulation order across channels can be randomized or not
nChannel = find(nOfm ~= 0);
if par.brand
    Order=nChannel(randperm(length(nChannel)));
else
    Order=nChannel; % Base to Apex, Apex to base, staggered.
end
StartPoint=vTimeOld(end)+Distance;
vTimeSlots=linspace(StartPoint,PulsePeriod+StartPoint,par.n); % Time vector only needed for plotting

%% Preallocation
ElectrodeScheme=zeros(par.m,par.n);
if any(nOfm) %If in the current time-slice there is a stimulation somewhere, then the n-of-M Order is important
for pc=1:par.n
        ElectrodeScheme(Order(pc),pc)=nOfm(Order(pc)); % Zuordnung der Amplitude auf Elektrodenplatz
end
end

% %Eliminate rows that have zeros.
% % only relevent when par.n < 12
% if any(ElectrodeScheme(:)) > 0 % are there amplitudes in frame > 0, else it is silence
%     ElectrodeScheme(:,~any(ElectrodeScheme,1))=[]; % eliminate rows that only contain zeros so only n channels are stimulated
    nSlots=par.n;
%     if nSlots~=par.n% this case should not happen
%         display('Error in GeneratePulseCycles')
%         keyboard
%     end
    [StimulationOrder, ~]=find(ElectrodeScheme>=0); %returns row = Order of stimulated electrodes
%     if isempty(StimulationOrder)
%         StimulationOrder = [1:par.n]; %It will only enter, if our signal has only silence in there. Then just assume one electrode after the other with zero amplitude each
%     end
% else % if there are only zeros in the signal frame
%     nSlots=0;
%     StimulationOrder=double.empty(1,0);
%     ElectrodeScheme=double.empty(1,0);
%     vTimeSlots = StartPoint;
% end






















