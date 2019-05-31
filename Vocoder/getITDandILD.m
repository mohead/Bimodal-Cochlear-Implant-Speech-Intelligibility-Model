function [itd, ild, freqs] = getITDandILD(signal,fs)
    % load auditory modeling toolbox
%     ltfatstart; amtstart;
   % run mathias IPD model:
   [hairc_fine, fc, hairc_ild, env]=dietz2011(signal,fs);
    
   % ectraxt ILD and ITD for this signal:
   itd = median(hairc_fine.itd);
   ild = median(hairc_ild);
   freqs = fc;
end