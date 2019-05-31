function [vFrame, vTimeFrame, FrameIdx,DiffLength] = myWindowing (vSignal,SamplingRate, FrameLength,FrameIdx)
%% function to conduct Block Processing
% INPUT
% ---------
% vSignal :      Signal vector
% SamplingRate : Sampling rate of signal
% FrameLength :  Length of each frame [s]
% FrameShift :   Shift of consecutive frames [s]
% FrameIdx :     First sample of current frame
%
% OUTPUT
% ---------
% mFrames:       Matrix with values of each frame
% vTimeFrame:    vector of times in center of each frame [s]
% NumOfFrames:   Number of Frames in mFrames
% DiffLength:    Length that signal is shorter than frame length

DiffLength=0;
FrameLength=FrameLength*SamplingRate; % transform FrameLength from s to samples

if length(vSignal) < FrameIdx+FrameLength % if size of signal not sufficient for framesize
    DiffLength=round((FrameIdx+FrameLength)-length(vSignal)); 
    vSignal=[vSignal; zeros(DiffLength,1)]; % Zero-Padd last frame   
end

EndFrame=(FrameIdx+FrameLength)-1; % last sample of frame

vFrame=vSignal(FrameIdx:EndFrame); % current Frame of Signal
vTimeFrame=(FrameLength/2+(FrameIdx-1))/SamplingRate; % center time of frame






