function [SRT,stData,frameData] = runBSIM(mfSignal, varargin)
%         Parameters:
%          'fs'
%          'audiogram'
%                 [frequency HL_left HL_right]
%          'data'
%          'level'
%          'frame'
%          'framesec'
%          'impfkt'
%          'errorflag'
%
% -------------------------------------------------------------------------
% constants (defaults, can be changed by call parameters)
% -------------------------------------------------------------------------

stControl     = ParParser(length(mfSignal),varargin{:});


%%
% It is assumed that the equalization parameters do not depend on the
% internal noise.
% Therefore, the internal noise is set to zero (hearing threshold = -Inf) for the calculation of the
% equalization parameters.(see stControl.HL.mfAudiogram in line 24)
% The hearing threshold is taken into consideration later
% in the Quick BSIM. Use the following mfAudiogram for your audiogram:
% mfAudiogram = [0 0 0; 50000 0 0]; % audiogram for SII and Quick-BSIM calculation
%%


% -------------------------------------------------------------------------
% processing starts here
% -------------------------------------------------------------------------
if ~stControl.run.bSkipBSIM
    
    % -------------------------------------------------------------------------
    % size check
    % -------------------------------------------------------------------------
    if size(mfSignal,2) ~= 4
        error('please provide 4 signal columns in the order of [speech_L speech_R noise_L noise_R]');
    end
    if size(mfSignal,1) < stControl.signal.iFrameLen
        error(['signal length is too short (must be > ' num2str(iFrameLen) ' samples @ 44100kHz)']);
    end
    
    vfLevel = max(10*log10(mean(abs(mfSignal).^2)));
    mfSignal = mfSignal .* 10.^(-vfLevel/20);
    iNumFrames = floor((size(mfSignal,1)-stControl.signal.iFrameLen)/stControl.signal.iFrameShift)+1;
    
    stControl.model.Delays = zeros(30,iNumFrames);
    stControl.model.Alphas = zeros(30,iNumFrames);
    
    stControl.model.Counter = 1;
    
    c = 0;
    tic;
    % generate random noise
    internal_noise = 1;
    if internal_noise
        InternalNoise = randn(length(mfSignal),2);
        InternalNoise(:,1) = InternalNoise(:,1)./rms(InternalNoise(:,1));
        InternalNoise(:,2) = InternalNoise(:,2)./rms(InternalNoise(:,2));
    else
        InternalNoise = zeros(length(mfSignal),2);
    end
    
    for iFrameIdx = 1:iNumFrames
        tmp = mfSignal((iFrameIdx-1)*stControl.signal.iFrameShift+(1:stControl.signal.iFrameLen),:).*repmat(stControl.signal.vWindow,1,4);
        IntNoi = InternalNoise((iFrameIdx-1)*stControl.signal.iFrameShift+(1:stControl.signal.iFrameLen),:).*repmat(stControl.signal.vWindow,1,2);
        % Start BSIM:
        [tmpRes,stControl] = BSIM03(tmp,IntNoi,stControl);
        %[tmp_Res,stControl] = BSIM_instant(mfSignal,stControl);
        fn = fieldnames(tmpRes);
        for k=1:length(fn);
            stResult.(fn{k})(:,:,iFrameIdx)=tmpRes.(fn{k});
        end
        
        if iNumFrames > 1 && iFrameIdx/iNumFrames > c
            c = c+0.2;
            elapsed = toc;
            rem_time = elapsed*(iNumFrames/iFrameIdx-1);
            disp(sprintf('estimated remaining time (%3.0f%% done): %02.0f:%02.0f',iFrameIdx/iNumFrames*100,floor(rem_time/60),mod(rem_time,60)));
        end
    end
end

SRT = NaN*ones(iNumFrames,1);
for iFrameIdx = 1:iNumFrames
    tmpRes.vfInt   = stResult.vfInt(:,:,iFrameIdx);
    tmpRes.vfCross = stResult.vfCross(:,:,iFrameIdx);
    tmpRes.Alpha   = stResult.Alpha(:,:,iFrameIdx);
    
    % Start Quick BSIM:
    [Sp,Ns,BetterEar]       = QuickBSIM(tmpRes,stControl,vfLevel*ones(1,4),stControl.run.audiogram_QCKBSIM );
    
    midslopevalue=stControl.Bimodal_SII_Switch_value;midslopevalue_inv=30-midslopevalue;
    a = sum(BetterEar==2);
    SII_act_func=(((a-midslopevalue)./(1+abs(a-midslopevalue)).*midslopevalue*(midslopevalue_inv+1))...
+(midslopevalue-1)*(midslopevalue_inv+1))...
/(midslopevalue*midslopevalue_inv+(midslopevalue-1)*(midslopevalue_inv+1));
    SII2IntellParam(:,2)    = stControl.run.SIIval(:,1).*(1-SII_act_func) + stControl.run.SIIval(:,2).*(SII_act_func); 
%     SII2IntellParam(:,2)    = (SIITemp(:,1).*sum(BetterEar==1)+SIITemp(:,2).*sum(BetterEar==2))./length(BetterEar);
    % Computation of Speech Reception Threshold (SRT) using SII
    for ll = 1:size(SII2IntellParam,1)
%         stResult.SRT(iFrameIdx,ll) = SRTFromSII(Sp(:),Ns(:),-Inf*ones(size(Ns(:))),nImpFkt,...
%             [stControl.model.FB.vfCenterFreqs(:), stControl.model.FB.vfBandWidth(:)],...
%             0,SII2IntellParam(ll,1),SII2IntellParam(ll,2));
        stResult.SRT(iFrameIdx,ll) = SRTFromSII(Sp(:),Ns(:),-Inf*ones(size(Ns(:))),stControl.model.ImpFkt,...
            zeros(30,1),... %insertion gain
            0,SII2IntellParam(ll,1),SII2IntellParam(ll,2));
        
        frameData.Sp(:,iFrameIdx)  = Sp(:);
        frameData.Ns(:,iFrameIdx)  = Ns(:);
    end
    
end

stData.stControl    = stControl;
stData.stResult     = stResult;
stData.iNumFrames   = iNumFrames;
stData.vfLevel      = vfLevel;

% % nanmean requires Statistics Toolbox
% SRT = nanmean(stResult.SRT);
% overall SRT is average over all filter
SRT = mean(stResult.SRT,1);
frameData.SRT = stResult.SRT;
frameData.vfLevel = vfLevel;
frameData.mfAudiogram = stControl.run.audiogram_QCKBSIM;
% end of file
end
 
