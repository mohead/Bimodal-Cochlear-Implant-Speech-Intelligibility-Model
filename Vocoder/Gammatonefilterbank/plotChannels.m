% Function to help with plotting differenct time signals in different channels from a channel x time
% matrix. It returns the rows and columns speperately and you can just
% plot() them.
function plotChannels(Data, fs, scalefactor, varargin)

[rows, columns] = getplotData(Data, fs, scalefactor);
plot(rows, columns, varargin{1:end-1});
title(varargin{end});

% Subfunction to calculate rows and columns
function [ax val] = getplotData(Data, fs, scalefactor)
[nRows nCols]=size(Data);
colValues=0:1/fs:(nCols-1)/fs; % Time axis
val = bsxfun(@plus,Data.*scalefactor,[1:nRows]');
% peakGain=nRows*waveHeight;
% a=max(max((Data)))*(0:nRows-1)';
% 
% % lower frequency channels mask higher frequency channels
% x=peakGain*Data+repmat(a,1,nCols);
% x=nRows*x/max(x, [], 2);
% 
% % for row=2:nRows
% %     x(row,:)=max(x(row-1,:), x(row,:));
% % end
ax = colValues;
% val = x';


end

end
% end of file