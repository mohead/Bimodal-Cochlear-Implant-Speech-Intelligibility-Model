% function [scaled_signal,gain_per_channel] = setSignalRMS(signaltoscale,RMSvalue,dBswitch)
% This function sets a signal to a given rms-value and returns the scaled signal.
% The rms-value can be either a 'linear' factor or in 'dB'.
% Specify this by setting the dbswitch accordingly.
% Additionally you can also return the gain per channel, to know how much a signal was scaled.
% The unit of the gain per channel is the same as specified in dBswitch.
%
% Example usages:
% Set signal to same rms as another signal:
% [scaled_signal] = setSignaltoRMS(signal_a,rms(signal_b),'linear');
% 
% Set signal to -30 dB FS:
% [scaled_signal, dB_gain_per_channel] = setSignaltoRMS(signal_a,-30,'dB');
% 
% This function relys on the rms-script written for psylab by Martin
% Hansen.
% For unit-Test look at Test_setSignalRMS.m
%
% Author: Ben Williges, ben.williges@uni-oldenburg.de
function [scaled_signal,gain_per_channel] = setSignaltoRMS(signaltoscale,RMSvalue,dBswitch)
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
       if size(signaltoscale,1) < size(signaltoscale,2)
           rms_current_signal = rms2(signaltoscale,2);
       else
           rms_current_signal = rms2(signaltoscale,1);
       end
       linear_gain_per_channel = bsxfun(@rdivide,RMSfactor, (rms_current_signal+eps));
       scaled_signal = bsxfun(@times,signaltoscale,linear_gain_per_channel); %Works also for matrixes
       
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