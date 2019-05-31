% % function [out_signal,mspatial_weight] = CIinteract_Ben(signal,m, distance_electrodes, lambda)
% The spatial spread is simulated by an exponentially decaying function:
% Amplitude.*exp^(-x/lambda), where x is the distance in mm from the stimulating electrode, lambda is the decaying constant of the electric field inside of the perylimph and the Amplitude is the current at the stimulating electrode.
% In the implementation the spatial spread is controlled by three variables:
% - lambda: This is the decaying constant in dB/m. You can find values for these in the litarature, e.g. Bingabr et al, 2008.
% - distance_electrodes: The distances in m between adjunct electrodes. E.g. cochlear has smaller distances between the electrodes than Med-El. This would be electrode-array dependent.
% - m: This influences, after how many electrodes you still want to calculate weights for the spatial spread (This depends on the total number of electrodes, choose them in a way, that you calculate the spatial spread along all electrodes, when you choose the electrode in the middle of the array).

function [out_signal,mspatial_weight] = CIinteract_Ben(signal,m, distance_electrodes, lambda)
% baue gewichtungsmatrix mit den Gewichtungen pro Kanal auf

d = distance_electrodes; %Distance between electrodes in [m]

% The spatial spread in dB/m is assumed to be lambda
% 
Nr_chans = size(signal,1);
% 
mspatial_weight = eye(Nr_chans,Nr_chans)./2; %Wir addieren nachher nochmals die matrix dazu, so dass die Gewichte stimmen

weights_right = exp(-[1:m]*d/lambda);
weights_left = fliplr(weights_right);

tic
for ii = 1:size(mspatial_weight,1)
    mspatial_weight(ii,ii+1:m+ii) = weights_right;
end
mspatial_weight = mspatial_weight(1:Nr_chans,1:Nr_chans);
% Spiegele oberen Teil der Matrix an Hauptdiagonaler in den unteren Teil
mspatial_weight=triu(mspatial_weight,-1)'+mspatial_weight; 

puls_idx = signal>0;
        out_signal = zeros(size(signal));
for ii = 1:size(signal,1); %Over all channels
    channel_active_idx = puls_idx(ii,:) > 0; % Finde alle Zeitpunkte, in denen der aktuelle Kanal aktiv ist
    out_signal(:,channel_active_idx) = bsxfun(@times, signal(ii,channel_active_idx),mspatial_weight(ii,:)');
end
end

