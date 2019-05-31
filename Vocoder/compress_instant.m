function outputsignal = compress_instant(inputsignal)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instantaneous compression, i.e. sample-by-sample of the amplitude of
% signal in different bands
% input: inputsignal with n bands across time t (n,t)
% output: outputsignal compressed with n bands across time t (n,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%These values are taken from the ACE-demo file.
B = 0.0156; %base level
M = 0.5859; %saturation level
alpha_c = 415.96; % controls the steepness of the compression function
load('lookupcompress.mat');

%initialize outputsignal
outputsignal = zeros(size(inputsignal));

%loop over frequency channels
for iCounter = 1:size(inputsignal,1)
    % this is the slow variant, which uses the log
    %outputsignal(iCounter,:) = process_compression_ci(inputsignal(iCounter,:),B,M,alpha_c);
    
    % this is the fast variant, which uses a lookup-table
    idx = round(1000*(inputsignal(iCounter,:))+1);
    outputsignal(iCounter,:) = out(idx);
end

