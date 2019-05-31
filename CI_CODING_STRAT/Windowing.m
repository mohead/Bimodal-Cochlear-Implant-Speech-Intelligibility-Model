function  [vFrame, vTimeFrame, FrameIdx,DiffLength]= Windowing (vSignal,SamplingRate, FrameLength,FrameIdx)
%% function to conduct Block Processing
% INPUT
% ---------
% vSignal :      Signal vector
% SamplingRate : Sampling rate of signal
% FrameLength :  Length of each frame [s]
% FrameShift :   Shift of consecutive frames [s]
%                
%
% OUTPUT
% ---------  
% mFrames:       Matrix with values of each frame
% vTimeFrame:    vector of times in center of each frame [s]
% NumOfFrames:   Number of Frames in mFrames
%
% Author:
% -----------
% Anja Eichenauer 



DiffLength=0;

if length(vSignal) < FrameIdx+FrameLength % Does last frame exceed signal length?
    DiffLength=round((FrameIdx+FrameLength)-length(vSignal)); 
    vSignal=[vSignal; zeros(DiffLength,1)]; % Zero-Padd last frame   
end

EndFrame=(FrameIdx+FrameLength)-1;

vFrame=vSignal(FrameIdx:EndFrame); % current Frame of Signal
vTimeFrame=(FrameLength/2+(FrameIdx-1))/SamplingRate; % center time of frame


