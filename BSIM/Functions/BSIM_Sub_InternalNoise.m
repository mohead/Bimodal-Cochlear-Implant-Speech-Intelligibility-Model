%% BSIM_Sub_InternalNoise --------------------------------------------
function [stControl, vfInternalNoise] = BSIM_Sub_InternalNoise(stControl)
%---------------------------------------------------------------
% function: BSIM_Sub_InternalNoise 
% calculates the internal noise intensities for the left and right ear according
% to the audiogram
%
%---------------------------------------------------------------

mfAudiogram     = stControl.HL.mfAudiogram;     % the individual audiogram
mfStandardAudio = stControl.HL.mfFreqHL2SPL;    % convert dB HL to dB SPL (internal noise spectrum level NH)

fFc_Hz          = stControl.model.FB.vfCenterFreqs(stControl.band.iCfIdx);

% interp1 cannot extrapolate linearly. use next valid HL instead, if center
% frequency is out of range (both audiogram and HL to SPL conversion)
fFc_Hz(fFc_Hz < max(min(mfAudiogram(:,1)),min(mfStandardAudio(:,1)))) = max([fFc_Hz(:); min(mfAudiogram(:,1)); min(mfStandardAudio(:,1))]);
fFc_Hz(fFc_Hz > min(max(mfAudiogram(:,1)),max(mfStandardAudio(:,1)))) = min([fFc_Hz(:); max(mfAudiogram(:,1)); max(mfStandardAudio(:,1))]);

% folgender Kommentarblock kann auch raus (Christopher):

% vfInternalNoise = stControl.signal.iSigLen/stControl.signal.iFs * ...
%     10.^(( ...
%          interp1(    mfAudiogram(:,1),    mfAudiogram(:,[2 3]),fFc_Hz,'linear','extrap')...
%         +interp1(mfStandardAudio(:,1),mfStandardAudio(:,[2 2]),fFc_Hz,'linear','extrap')...
%     )/10)/(10^(4/10)-1);

% changed reference to match approximately SII reference internal noise
% spectrum, see also interneal noise calc procedure

vfInternalNoise = stControl.signal.iSigLen/stControl.signal.iFs...
                .*repmat(stControl.model.FB.vfBandWidth(stControl.band.iCfIdx),1,2)...
                .* ...
    10.^((stControl.HL.fdBDetThresh + ...
         interp1(    mfAudiogram(:,1),    mfAudiogram(:,[2 3]),fFc_Hz,'linear','extrap')...
        +interp1(mfStandardAudio(:,1),mfStandardAudio(:,[2 2]),fFc_Hz,'linear','extrap')...
    )/10);

vfInternalNoise(isnan(vfInternalNoise)) = 0; % if NAN, set internal noise to zero
% end of file