function [Speech,Noise,BestEar] = QuickBSIM(stResult,stControl,vfNoiseLevels,mfAudiogram)

stControl.HL.mfAudiogram = mfAudiogram;
stControl.band.iCfIdx = 1:length(stControl.model.FB.vfCenterFreqs);

% calculate Internal noise (taking the audiogram into account)
[stControl,HTL] = BSIM_Sub_InternalNoise(stControl);

gain = 10.^(stResult.Alpha*[1 -1]/10);
gain(gain > 1) = 1;

vfLevels = 10.^(vfNoiseLevels/10);

IS(:,1) = vfLevels(1).*stResult.vfInt(:,1);
IS(:,2) = vfLevels(2).*stResult.vfInt(:,2);
IN(:,1) = vfLevels(3).*stResult.vfInt(:,3) + HTL(:,1); % add internal noise to external noise (left ear)
IN(:,2) = vfLevels(4).*stResult.vfInt(:,4) + HTL(:,2); % add internal noise to external noise (right ear)

if stControl.model.EC.EC_flank==1
    S =   gain(:,1).*IS(:,1) ...
        + gain(:,2).*IS(:,2) ...
        + sqrt(prod(gain,2)*prod(vfLevels([1 2]))).* stResult.vfCross(:,1);
    N =   gain(:,1).*IN(:,1)...
        + gain(:,2).*IN(:,2)...
        + sqrt(prod(gain,2)*prod(vfLevels([3 4]))).* stResult.vfCross(:,2);
    % Calculate monaural SNR
    Mono = [IS./IN];
    % find better ear
    [MaxMono,Ch]    = max(Mono,[],2);
    BestEar         = Ch;             % indicate the best ear (use best ear for binaural too)
    % find idx where binaural SNR is below the monaural SNR (better ear listening)
    idx = find(S./N < MaxMono); %this is the switch between SNR-monaural and SNR-binaural!!
    S(idx) = IS(sub2ind(size(IS),idx,Ch(idx)));
    N(idx) = IN(sub2ind(size(IN),idx,Ch(idx)));
    % Calculate Speech and Noise Level
    Speech = (10*log10(S*stControl.signal.iFs/stControl.signal.iSigLen./stControl.model.FB.vfBandWidth));
    Noise  = 10*log10(N*stControl.signal.iFs/stControl.signal.iSigLen./stControl.model.FB.vfBandWidth);
else
    
    if stControl.model.EC.BetterEar==1
        % Calculate monaural SNR
        Mono = [IS./IN];
        % find better ear
        [MaxMono,Ch] = max(Mono,[],2);
        BestEar         = Ch;     	% indicate the best ear (use best ear for binaural too)
        % find idx where binaural SNR is below the monaural SNR (better ear
        % listening)
        for idx =1:30;
            S(idx,1) = IS(idx,Ch(idx));
            N(idx,1) = IN(idx,Ch(idx));
        end
     % Calculate Speech and Noise Level
        Speech = (10*log10(S*stControl.signal.iFs/stControl.signal.iSigLen./stControl.model.FB.vfBandWidth));
        Noise  = 10*log10(N*stControl.signal.iFs/stControl.signal.iSigLen./stControl.model.FB.vfBandWidth);
    
    elseif stControl.model.EC.BetterEar==2
        S = IS(:,1);
        N = IN(:,1);
        Speech = (10*log10(S*stControl.signal.iFs/stControl.signal.iSigLen./stControl.model.FB.vfBandWidth));
        Noise  = 10*log10(N*stControl.signal.iFs/stControl.signal.iSigLen./stControl.model.FB.vfBandWidth);
        BestEar = ones(30,1);    % Left Ear always used 
    
    elseif stControl.model.EC.BetterEar==3
        S = IS(:,2);
        N = IN(:,2);
        Speech = (10*log10(S*stControl.signal.iFs/stControl.signal.iSigLen./stControl.model.FB.vfBandWidth));
        Noise  = 10*log10(N*stControl.signal.iFs/stControl.signal.iSigLen./stControl.model.FB.vfBandWidth);
        BestEar = ones(30,1).*2; % Right Ear always used 
    end
end
end
% end of file