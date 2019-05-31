function [vSpec vFreqSpec] = CISpectrum(vSignal,par)
%% function to perform FFT on signal
% INPUT
% ---------
% vSignal :      Signal vector
% par     :      Externally defined parameters
%
% OUTPUT
% ---------
% vSpec     :   frequency magnitudes (Spectrum)
% vFreqSpec :   frequency vector [Hz]

fs=par.CIFs;
NFFT=128;
vSpec=fft(vSignal,NFFT);
LSpec=length(vSpec);
vSpec=vSpec(1:LSpec/2+1); % only keep half the spectrum + Nyquist frequency
Freq=0:LSpec/2;
vFreqSpec= Freq./length(vSignal).*fs; % frequency vector