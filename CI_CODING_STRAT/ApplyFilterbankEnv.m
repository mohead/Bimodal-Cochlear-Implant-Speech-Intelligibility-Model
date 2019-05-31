function mSigCh = ApplyFilterbankEnv(bins, vSpec)
%% function to apply 22 channel FFT filter bank on signal and calc envelope per channel
% INPUT
% ---------
% bins     :     Matrix of zeros and ones. Ones at bins of channels that
%                are summarized. Size 22x65 if NFFT = 128.
% vSpec    :     Spectrum of acoustic signal
%
% OUTPUT
% ---------
% mSigCh    :     Matrix with weighted envelope magnitude per channel


Weights=[0 0.98 0.68 0.65];
vWeights=[Weights(1)*ones(2,1); Weights(2)*ones(9*1,1);...
    Weights(3)*ones(4*2,1); Weights(4)*ones(2*3+2*4+2*5+6+7+8,1); 0]; % Weight per frequency bin. The numbers in the ones are multiplied by
% the fft-bins inside of the corresponding channels. Last bin is set to
% zero, because it contains the DC-Part, which we don't want to stimulate

% Power-sum combination of bins (see Swanson formular 4.2.3)
mSigCh=sqrt(bins*(vWeights.*abs(vSpec).^2)); % Matrixmultiplication sums all amplitudes of Spec with respect to Ch
