function delay = Gfb_Delay_new(analyzer, delay_samples)
% delay = Gfb_Delay_new(analyzer, delay_samples)
%
% Gfb_Delay_new creates a new Gfb_Delay object that can act as the first stage
% of a synthsizer that resynthesizes the output of the gammatone filterbank
% analyzer.  The purpose of the delay object is to delay the output of each
% channel by a channel-dependent ammount of samples, so that the envelope of
% the impulse response of the analyzer is as large as possible at the desired
% delay.
% Additionally, the delay object will multiply this delayed output with a
% channel-dependent complex phase factor, so that the real part of the impulse
% response has a local maximum at the desired delay.  Finally, the delay ob-
% ject will output only the real part of each channel.
%
% The phase factors are approximated numerically in this constructor.  The
% approximation assumes parabolic behaviour of the real part of the impulse
% response in the region of the desired local maximum:  The phase factors
% are chosen so that the real parts of the impulse response in the samples
% directly preceeding and following the desired local maximum will be equal
% after multiplication with the pase factor.
%
% PARAMETERS:
% analyzer      The Gfb_Analyzer structure as returned by Gfb_Analyzer_new.
%               May also be a structure array of several pluggable analyzers.
% delay_samples The desired group delay in samples. must be at least 1,
%               because of the way the phase factors are computed.  Larger
%               delays lead to better signal quality
% delay         The new Gfb_Delay object
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002 Nov 2003

% filename : Gfb_Delay_new.m


delay.type           = 'Gfb_Delay';
analyzer             = Gfb_Analyzer_clear_state(analyzer);
impulse              = zeros(1, delay_samples + 2);
impulse(1)           = 1;
number_of_channels   = length(analyzer(1).center_frequencies_hz);

impulse_response     = impulse;
for a = analyzer
    impulse_response = Gfb_Analyzer_process(a, impulse_response);
end
[dummy, max_indices] = max(abs(impulse_response).');
max_indices          = max_indices - (max_indices > delay_samples + 1);

delay.delays_samples = delay_samples + 1 - max_indices;

delay.memory         = zeros(number_of_channels, ...
                             max(delay.delays_samples));

slopes = zeros(1, number_of_channels);
for channel = [1:number_of_channels]
  channel_max_index = max_indices(channel);
  slopes(channel) = (impulse_response(channel, channel_max_index+1) - ...
                     impulse_response(channel, channel_max_index-1));
end
slopes = slopes ./ abs(slopes);
delay.phase_factors = 1i ./ slopes;


%%-----------------------------------------------------------------------------
%%
%%   Copyright (C) 2002 2003  AG Medizinische Physik,
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
