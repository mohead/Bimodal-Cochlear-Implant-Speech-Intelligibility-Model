%% Function to calculate the timing differences between sequentially electrodogramm
% pulses. It takes as input an channel x timeframes matrix and the sampling
% frequency of the timeframes.
% As an ouput it returns a squeezed electrodogramm (e.g. timeframes, where 
% every channel is zero, gets removed), an vector of the stimulating
% electrodes (same time-dimension as the squeezed electrodogramm), an the
% corresponding time vector for the stimulation of electrodes.
% Caution: This function assumes, that there is ONLY sequential
% stimulation, e.g. all simultaneous stimulation is removed beforehand.
% Also it assumes, that there are only positive pulses, each of which is
% long enough to include positive phase, negative phase and inter-pulse-gap
% between the positive and negative phases.
function [squeezed_electrodogramm, vElectrodes, vDistance] = calculateDistance(electrodogramm, fs, pulselength,ipg)
    % go through electrodogramm sample by sample:
    vDistance = 0;
    vElectrodes = 0;
    squeezed_electrodogramm = [];
    [channels, sample_idx] = ind2sub(size(electrodogramm),find(electrodogramm>0));
    for ii = 1:length(sample_idx)
            dist = sample_idx(ii)/fs-sum(vDistance); % distance in seconds
            
                if vElectrodes(end) == channels(ii) && electrodogramm(channels(ii),sample_idx(ii)-1)>0
                % do nothing, if we operate on the tail of the pulse (e.g. pulselength >1 sample 
                else
                    if dist > 2*pulselength+ipg+4e-6 % pulslength check
                    vDistance = [vDistance dist];
                    vElectrodes = [vElectrodes channels(ii)];
                    squeezed_electrodogramm = [squeezed_electrodogramm electrodogramm(:,sample_idx(ii))];
                    end  
                end
    end
    vDistance = vDistance(2:end); % Remove zero from beginning
    vElectrodes = vElectrodes(2:end); %Remove zero from beginning
    assert(min(vDistance) > 2*pulselength+ipg+4e-6, 'Timing calculation is wrong');
end