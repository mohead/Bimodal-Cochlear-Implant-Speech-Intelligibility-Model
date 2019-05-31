% function [scaled_signal] = setSignalRMS(signaltoscale,RMSvalue,dBswitch)
% This function scales a signal by a given rms-value and returns the scaled signal.
% The rms-value can be either a 'linear' factor or in 'dB'.
% Specify this by setting the dbswitch accordingly.
%
% Example usages:
% Scale signal by 3 dB:
% [scaled_signal] = scaleSignalbyRMS(signal_a,3,'dB');
% 
% Scale signal by -sqare root(2):
% [scaled_signal] = scaleSignalbyRMS(signal_a,-sqrt(2),'linear');
% 
% This function relys on the rms-script written for psylab by Martin
% Hansen.
% For unit-Test look at Test_scaleSignalbyRMS.m
%
% Author: Ben Williges, ben.williges@uni-oldenburg.de
function [scaled_signal] = scaleSignalbyRMS(signaltoscale,RMSvalue,dBswitch)
       %% Input check
       switch dBswitch
           case 'dB'
               RMSfactor = 10.^(RMSvalue ./ 20);
           case 'linear'
               RMSfactor = RMSvalue;
           otherwise
               error('You must specify wether the RMS Value is ''linear'' or ''dB''');
       end
       
       %% Scaling
       scaled_signal = bsxfun(@times,signaltoscale,RMSfactor); %Works also for matrixes

       %% output arguments checks
       if nargout > 1
           switch dBswitch 
               case 'dB'
               gain_per_channel = 20 .* log10(linear_gain_per_channel+eps);
               case 'linear'
               gain_per_channel = linear_gain_per_channel;
               otherwise
               error('You must specify wether the RMS Value is ''linear'' or ''dB''');    
           end
       end
end