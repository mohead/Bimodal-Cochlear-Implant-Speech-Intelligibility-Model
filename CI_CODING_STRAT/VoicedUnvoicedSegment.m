function [VoicedUnvoiced, NZero] = VoicedUnvoicedSegment(vSignal)
%% function to decide if Voiced or Unvoiced via Zero crossings
%% returns 1 for voiced segments and 0 for unvoiced
% INPUT
% ---------
% vSignal:       One Frame of Signal
%
%
% OUTPUT
% ---------
% VoicedUnvoiced: Voiced or Unvoiced (1 for voiced segments and 0 for unvoiced)
% NZero:          Zero Crossing rate
%
% Author:
% -----------
% Anja Eichenauer

% Energy Detection
if sqrt(mean(vSignal.^2)) < 0.0156 % 0.0156 == par.B
    VoicedUnvoiced=0; % not enough energy in signal = unvoiced
    NZero=0;
else
    ZCLimit=0.2; % normalized Limit for estimation voiced/unvoiced
    posAmp=vSignal>0; % find idx with amplitudes > 0 (result: 0 or 1) %% Hier auf Breitbandsignal, bei FSP schmalband
    FindCrossings=diff(posAmp); % differences not equal to 0 == crossings
    NCrossings = length(find(FindCrossings)); % counts all nonzero elements == crossings
    NZero=NCrossings/length(vSignal); % normalized Number of Zero-Crossings
    
%     decide Voiced / unvoiced
    if NZero > ZCLimit
        VoicedUnvoiced=0; % unvoiced
    else
        VoicedUnvoiced=1; % voiced
    end
    
end



