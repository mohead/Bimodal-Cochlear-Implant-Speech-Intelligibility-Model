function f0= fZeroEstimator(Fs, vFrame)
%% function to estimate fundamental frequency
% INPUT
% ---------
% Fs :           Sampling rate of signal
% vFrame:        signal vector 
%                
%
% OUTPUT
% ---------  
% f0:            fundamental frequency
%
% Author:
% -----------
% Anja Eichenauer 


rxx=xcorr(vFrame); % autocorrelation of frame
StartAnalysis=round((length(rxx)+1)/2); % whole length of conv is N+M-1, we only need half of it
% Only search for f0 within a limited time frame
EstimationAreaStart=250; %Hz (limit maximal f0)
EstimationAreaEnd=80; %Hz (limit lowest f0)
AnalysisWindowEnd=round((1/EstimationAreaEnd)*Fs); % Time at 80 Hz
AnalysisWindowStart=round((1/EstimationAreaStart)*Fs); % Time at 250 Hz

Start=StartAnalysis+AnalysisWindowStart; % Thats the Area the analysis actually starts
Ende=StartAnalysis+AnalysisWindowEnd; % And where it ends
rxxAnalysis=rxx(Start:Ende); % ACF 

% find Maximum within analysis window (which equals T@f0)
[c i]=max(rxxAnalysis);
Trxx=((i-1)+AnalysisWindowStart)/Fs; % compute T from max, subtract 1 because it doesn´t start at t=0s
f0=1/Trxx; % compute f0 from T