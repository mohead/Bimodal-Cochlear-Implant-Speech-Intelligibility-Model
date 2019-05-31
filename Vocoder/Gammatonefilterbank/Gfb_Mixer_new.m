function mixer = Gfb_Mixer_new(analyzer, delay, samples, iterations)
% mixer = Gfb_Mixer_new(analyzer, delay, samples, iterations)
% 
% Gfb_Mixer_new creates a Gfb_Mixer object with gain factors suitable
% to calculate a weighted sum of the channels present in the output of the
% given delay.  The gain factors are computed numerically from an impulse
% response of <samples> duration.
% The <samples> and <iterations> arguments may be omitted, in which case
% reasonable default arguments are assumed.
%
% PARAMETERS
% analyzer   A Gfb_Analyzer structure as created by Gfb_Analyzer_new. The
%            mixer created by this function can act as part of a synthesizer
%            that resynthesizes the output of this analyzer
% delay      A Gfb_Delay structure as created by Gfb_Delay_new, Together with
%            the mixer created by this function, this delay can form a
%            synthesizer that resynthesizes the output of the analyzer
% samples    The number of samples of the impulse response that will be
%            used for calculating the spectrum of the resynthesized signal.
%            If this parameter is omitted, then samples will be chosen so
%            that the impulse response will be long enough to hold at least
%            GFB_GAINCALC_MIN_PERIODS (see Gfb_set_constants.m, usually =20)
%            full periods of the lower cutoff frequency of the analyzer
% iterations The gain factors are approximated numerically in iterations.
%            If this parameter is omitted, then the number of iterations will
%            be  GFB_GAINCALC_ITERATIONS (see Gfb_set_constants.m, usually
%            =100)
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002 Nov 2003
% Modified for Braecker-Vocoder: BW, Jan, 2015

% filename : Gfb_Mixer_new.m


global GFB_GAINCALC_MIN_PERIODS GFB_GAINCALC_ITERATIONS;
Gfb_set_constants;

mixer.type           = 'Gfb_Mixer';
analyzer             = Gfb_Analyzer_clear_state(analyzer);
delay                = Gfb_Delay_clear_state(delay);
center_frequencies   = analyzer(1).center_frequencies_hz;
number_of_channels   = length(center_frequencies);
sampling_frequency   = analyzer(1).sampling_frequency_hz;
order                = analyzer(1).gamma_order;
bw_factor            = analyzer(1).bandwidth_factor;
if length(analyzer) > 1
  if ~isequal(analyzer.center_frequencies_hz)
    error('All analyzers must have the same center frequencies');
  end
  if ~isequal(analyzer.sampling_frequency_hz)
    error('All analyzers must have the same sampling frequency');
  end
end

if (nargin < 3) % samples not given
    %Calculate lowest cutoff frequency
    % Calculate ERB (eq. 13)
    ERB = GFB_L + analyzer(1).center_frequencies_hz(1) / GFB_Q;
    % calculate bandwith (eq. 15 with term a_y from eq. 14)
    a_y = (pi * factorial(2*order-2) * 2^(-1*(2*order-2))) / ...
                   (factorial(order-1)^2);
    c_y = 2 * sqrt(2^(1 / order) -1);
    bandwith = c_y / a_y *ERB; % (eq. 15)
    lower_cutoff_frequency_hz = analyzer(1).center_frequencies_hz(1) - floor((bw_factor(1)*bandwith)/2);
  samples = 2^ceil(log2(GFB_GAINCALC_MIN_PERIODS * sampling_frequency / ...
			lower_cutoff_frequency_hz));
end
impulse              = zeros(1, samples);
impulse(1)           = 1;

% The indices needed to extract the values for the center_frequencies
% from an FFT vector.
spectrum_indices     = round(center_frequencies * samples ...
			     / sampling_frequency);
mixer.gains          = ones(number_of_channels, 1);

impulse_response     = impulse;
for a = analyzer
    impulse_response = Gfb_Analyzer_process(a, impulse_response);
end
impulse_response     = Gfb_Delay_process(delay, impulse_response);

response_spectrum    = fft(real(impulse_response).');
% 1st channel's spectrum is in column 1 of response_spectrum

if (nargin < 4)
  iterations = GFB_GAINCALC_ITERATIONS;
end
for i = [1:iterations]
  selected_spectrum = response_spectrum(spectrum_indices,:);
  % 1st channel's spectrum is still in column 1 of response_spectrum

  % add selected spectrum of all channels with gain factors
  selected_spectrum = selected_spectrum * mixer.gains;

  % calculate better gain factors from result
  mixer.gains = mixer.gains ./ abs(selected_spectrum);
end
mixer.gains = mixer.gains.';
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
