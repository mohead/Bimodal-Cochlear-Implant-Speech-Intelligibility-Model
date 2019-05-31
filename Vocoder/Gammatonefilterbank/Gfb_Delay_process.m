function [output, delay] = Gfb_Delay_process(delay, input)
% [output, delay] = Gfb_Delay_process(delay, input)
%
% Each channel (row) of the input data will be delayed by a channel-dependend
% ammount of samples, then multiplied with a channel-dependend complex
% constant.  Finally, the real part of this product will be returned.
%
% PARAMETERS
% delay   A Gfb_Delay structure created from Gfb_Delay_new.  The delay
%         will be returned with updated delayline states as the second
%         return parameter
% input   A complex matrix containing the signal to delay.  Each row
%         corresponds to a filterbank channel
% output  A real matrix containing the delay's output
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002

% filename : Gfb_Delay_process.m


[number_of_channels, number_of_samples] = size(input);
if (number_of_channels ~= length(delay.delays_samples))
  error('input rows must match the number of channels');
end
output = zeros(number_of_channels, number_of_samples);
for channel = [1:number_of_channels]
  if (delay.delays_samples(channel) == 0)
    output(channel,:) = ...
        real(input(channel,:) * delay.phase_factors(channel));
  else
    tmp_out = [delay.memory(channel,1:delay.delays_samples(channel)), ...
               real(input(channel,:) * delay.phase_factors(channel))];
    delay.memory(channel,1:delay.delays_samples(channel)) = ...
        tmp_out(number_of_samples+1:length(tmp_out));
    output(channel,:) = tmp_out(1:number_of_samples);
  end
end


%%-----------------------------------------------------------------------------
%%
%%   Copyright (C) 2002   AG Medizinische Physik,
%%                        Universitaet Oldenburg, Germany
%%                        http://www.physik.uni-oldenburg.de/docs/medi
%%
%%   Permission to use, copy, and distribute this software/file and its
%%   documentation for any purpose without permission by UNIVERSITAET OLDENBURG
%%   is not granted.
%%   
%%   Permission to use this software for academic purposes is generally
%%   granted.
%%
%%   Permission to modify the software is granted, but not the right to
%%   distribute the modified code.
%%
%%   This software is provided "as is" without expressed or implied warranty.
%%
%%   Author: Tobias Peters (tobias@medi.physik.uni-oldenburg.de)
%%
%%-----------------------------------------------------------------------------
