function [db]=pegel(Y,ref,window_length)
  if nargin < 3
    window_length = 0;
  end
 if nargin < 2
    ref =  1;
 end
 
 % remove windowed signal from level calculation
 Y = Y(1+window_length:end-window_length,:);

len=length(Y);
prms=sqrt((1/len)*sum(Y.^2));

if prms~=0
    db=20*log10(prms/ref);
else
    db=-120;
end
