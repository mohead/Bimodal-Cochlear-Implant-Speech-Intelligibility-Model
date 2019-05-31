function VerticalStep = PlotElectrodogram(mInput, vtime, title_string)
if nargin < 3
    title_string = '';
end

VerticalStep=1200;
figure;
VerticalPosition=0;
for iCount = 1:size(mInput,1)
    plot(vtime, mInput(iCount,:)+ VerticalPosition);
    VerticalPosition =  VerticalPosition + VerticalStep;
    hold on;
end

xlabel('Time (s)')
set(gca,'YTick',[0:VerticalStep:(size(mInput,1)-1)*VerticalStep]);
set(gca,'YTickLabel',{[1:1:size(mInput,1)]});
ylabel('Electrode number');
title(title_string);