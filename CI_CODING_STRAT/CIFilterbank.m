function [bins, CenterFreq] = CIFilterbank(par)
%% function to generate 12 channel FFT filter bank
% INPUT
% ---------
% par     :      Externally defined parameters
%
% OUTPUT
% ---------
% bins     :     Matrix of zeros and ones. Ones at bins of channels that
%                are summarized. Size 22x65 if NFFT = 128, e.g. (128/2 (symetrical spectrum) + 1 Bin DC-Part).

CenterFreq=[250 375 500 625 750 875 1000 1125 1250 1437 1686 1937 2187 2500 2875 3312 3812 4375 5000 5687 6500 7437]; % Cochlear Center Frequencies
NoOfCh=length(CenterFreq);
if NoOfCh ~= par.m
    error('Wrong filterbank configuration')
end
% How many bins fall into each channel:
NoOfBinsPerCh=[0 0 1 1 1 1 1 1 1 1 1 2 2 2 2 3 3 4 4 5 5 6 7 8]; % Total 62 bins
BinIdxAll=cumsum(NoOfBinsPerCh)+2;
BinIdx=BinIdxAll(3:end); % Brett Swanson PhD, p. 59:
%     In the SPrint and Freedom, M = 128, and the input sample rate is (approximately) 16 kHz. Bin 0 (fc = 0 Hz) and bin 64 (fc = 8000 Hz) 
% are purely real, and bins 1 to 63 are complex. Because the input signal is real, the output has Hermitian symmetry, 
% and bins 65 to 127 are not required. Thus there are 65 band-pass filters, with centre frequencies spaced linearly at 
% multiples of 125 Hz. For the Hann window, the -6 dB bandwidth is 250 Hz (2 bins) (Harris 1978).
% The FFT creates a bank of FIR filters with a linear frequency spacing of 125 Hz. 
% However, a linear-log frequency spacing is required (§4.3.1). For the bands below 1000 Hz, each band can be simply assigned to one FFT bin, 
% starting at bin 2 (centred at 250 Hz). Bins 0 and 1 are discarded. 
% For the higher frequency bands, two or more consecutive FFT bins are combined to produce 
% successively wider bands. There are two methods of combining bins: vector-sum and power-sum.

% generate bin matrix
bins=zeros(NoOfCh,BinIdxAll(end)+1); %+1 to account for the DC-Bin, which is the last one
StartBin=3;
for ii=1:NoOfCh
    bins(ii,StartBin:BinIdx(ii))= 1; 
    StartBin=BinIdx(ii)+1;
end
