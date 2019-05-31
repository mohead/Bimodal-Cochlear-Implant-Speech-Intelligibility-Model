 function plotaudi(vFreq, mDBHL, vSide, vTitle);
%
% function plotaudi(vFreq, mDBHL, vSide, vTitle);
%
% plot an audiogram.
% vFreq : frequency scale
% mDBHL : matrix of thresholds (each column corresponds to one side)
% vSide : vector describing which column in 'mDBHL' stands for which
%         side. EX.: 'rl' = first column is right, second is left.
%         'rl' is the default setting when 'mDBHL' is a 2 column matrix.
% vTitle: Title for the plot
%
% author/date: ja / 5.97
% modified by tb / 3.14
%

% set defaults for ungiven parameters
[z s] = size(mDBHL);
if nargin < 3,
	if s == 2,
		vSide = 'rl';
	else
		vSide = 'x';
	end;
end;
if nargin < 4,
	vTitle = 'Audiogram';
end;
	
% struct describing plot configuration
sTable.vFreq      = [125 250 500 750 1000 1500 2000 3000 4000 6000 8000]';
sTable.vTicks     = [1   2   3   3.5   4  4.5    5  5.5    6  6.5    7]';
sTable.vAxTicks   = [0.8 7.2 -10 130];
sTable.StyleNone  = 'y+-';
sTable.StyleRight = 'ro-';
sTable.StyleLeft  = 'bx-';

% struct desribing data
sData.vFreq     = vFreq(:);
sData.vDBHLLeft = [];
sData.vDBHLRight= [];
sData.vDBHLNone = [];
for i = 1:length(vSide),
	if 		(vSide(i) == 'l'),	sData.vDBHLLeft  = mDBHL(:,i);
	elseif	(vSide(i) == 'r'),	sData.vDBHLRight = mDBHL(:,i);
	else		sData.vDBHLNone = mDBHL(:,i);
	end;
end;
sData.vTicks    = interp1(sTable.vFreq, sTable.vTicks, sData.vFreq, 'linear');

% plot data and set axis
%fig;
hold on;
if (~isempty(sData.vDBHLLeft))
	plot(sData.vTicks, sData.vDBHLLeft, sTable.StyleLeft);
end;
if (~isempty(sData.vDBHLRight))
	plot(sData.vTicks, sData.vDBHLRight, sTable.StyleRight);
end;
if (~isempty(sData.vDBHLNone))
	plot(sData.vTicks, sData.vDBHLNone, sTable.StyleNone);
end;
axis(sTable.vAxTicks);
set(gca,'YDir','reverse');
set(gca,'XTick',sTable.vTicks);
set(gca,'XTickLabel',sTable.vFreq/1000);
set(gca,'XAxisLocation','top');
grid on;

% title and outfit
xlabel('Frequency [kHz]');
ylabel('Threshold [dBHL]');
h = line(sTable.vAxTicks(1:2),[0 0]);
set(h,'LineWidth',1);
set(h,'Color','black');



