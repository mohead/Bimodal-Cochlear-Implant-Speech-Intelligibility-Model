function analyzer = Gfb_Analyzer_new(sampling_frequency_hz,         ...
                                     lower_cutoff_frequency_hz,     ...
                                     specified_center_frequency_hz, ...
                                     upper_cutoff_frequency_hz,     ...
                                     filters_per_ERBaud,            ...
                                     gamma_order,                   ...
                                     bandwidth_factor,              ...
                                     normalization_divisor)
% analyzer = Gfb_Analyzer_new(sampling_frequency_hz,         ...
%                             lower_cutoff_frequency_hz,     ...
%                             specified_center_frequency_hz, ...
%                             upper_cutoff_frequency_hz,     ...
%                             filters_per_ERBaud             ...
%                             gamma_order,                   ...
%                             bandwidth_factor,              ...
%                             normalization_divisor)
%
% Gfb_Analyzer_new constructs a new Gfb_Analyzer object.  The analyzer
% implements the analysis part of a gammatone filterbank as described
% in [Hohmann 2002].
% It consists of several all-pole gammatone filters; each
% one with a bandwidth of 1 ERBaud (times bandwidth_factor),
% and an order of gamma_order.
% The center frequencies of the individual filters are computed as
% described in section 3 of [Hohmann 2002].
%
% PARAMETERS: (all frequencies in Hz)
% sampling_frequency_hz      The sampling frequency of the signals on which
%                            the analyzer will operate
% lower_cutoff_frequency_hz  The lowest possible center frequency of a
%                            contained gammatone filter
% specified_center_frequency_hz       ( == "base frequency")
%                            One of the gammatone filters of the analyzer
%                            will have this center frequency.  Must be >=
%                            lower_cutoff_frequency_hz
% upper_cutoff_frequency_hz  The highest possible center frequency of a
%                            contained gammatone filter.  Must be >=
%                            specified_center_frequency_hz
% filters_per_ERBaud         The density of gammatone filters on the ERB
%                            scale.
% gamma_order                optional:
%                            The order of the gammatone filters in this
%                            filterbank.
%                            If unspecified, the default value from
%                            Gfb_set_constants.m is used.
%                            A vector of multiple orders may be given.
%                            In this case, a structure array of multiple
%                            analyzers is returned.
% bandwidth_factor           optional:
%                            The bandwidth parameter of the individual filters
%                            is calculated from the Equivalent Rectangular
%                            Bandwidth (ERB) according to equation 14 in
%                            [Hohmann 2002]. ERB is taken from the Glasberg &
%                            Moore formula for a specific center frequency
%                            (equation 13 in [Hohmann 2002]).
%                            Using this parameter, it is possible to widen or
%                            narrow all filters of the filterbank with a
%                            constant bandwidth factor.
%                            A vector of multiple bandwidth factors may be
%                            given with the same size as gamma_order.  If 
%                            gamma_order is a vector and bandwidth_factor is a
%                            scalar, the same bandwidth_factor will be used for
%                            every created Gfb_Analyzer.
%                            The default value for this parameter is 1.
% normalization_divisor      optional parameter for internal correction of
%                            normalization_factor, do not use.
%
% OUTPUT:
% analyzer                   The constructed Gfb_Analyzer object(s).
%                            Multiple analyzers returned by this function
%                            in a structure array (when gamma_order is a
%                            vector) are pluggable when order of application
%                            is preserved. The first returned filter works for
%                            real input, while the remaining filters require
%                            the complex output from the previous filterbank.
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Jan, Sep 2003

% filename : Gfb_Analyzer_new.m

if (nargin < 6)
  % The order of the gammatone filter is derived from the global constant
  % GFB_PREFERED_GAMMA_ORDER defined in "Gfb_set_constants.m".  Usually,
  % this is equal to 4.
  global GFB_PREFERED_GAMMA_ORDER;
  Gfb_set_constants;
  gamma_order = GFB_PREFERED_GAMMA_ORDER;
end
if (nargin < 7)
  bandwidth_factor = 1;
end
if length(bandwidth_factor) == 1
  bandwidth_factor = bandwidth_factor * ones(size(gamma_order));
end
if (nargin < 8)
  normalization_divisor = 1;
end


if (length(gamma_order) ~= 1)
  for index = 1:length(gamma_order)
    analyzer(index) = Gfb_Analyzer_new(sampling_frequency_hz,         ...
                                       lower_cutoff_frequency_hz,     ...
                                       specified_center_frequency_hz, ...
                                       upper_cutoff_frequency_hz,     ...
                                       filters_per_ERBaud,            ...
                                       gamma_order(index),            ...
                                       bandwidth_factor(index),       ...
                                       normalization_divisor);
    normalization_divisor = 2;
  end
  return;
end

% To avoid storing information in global variables, we use Matlab
% structures:
analyzer.type                          = 'Gfb_Analyzer';
analyzer.sampling_frequency_hz         = sampling_frequency_hz;
analyzer.lower_cutoff_frequency_hz     = lower_cutoff_frequency_hz;
analyzer.specified_center_frequency_hz = specified_center_frequency_hz;
analyzer.upper_cutoff_frequency_hz     = upper_cutoff_frequency_hz;
analyzer.filters_per_ERBaud            = filters_per_ERBaud;
analyzer.fast                          = 0;


%
analyzer.center_frequencies_hz = ...
    Gfb_center_frequencies(filters_per_ERBaud, ...
			   lower_cutoff_frequency_hz,     ...
			   specified_center_frequency_hz, ...
			   upper_cutoff_frequency_hz);

% This loop actually creates the gammatone filters:
for channel = [1:length(analyzer.center_frequencies_hz)]
  center_frequency_hz = analyzer.center_frequencies_hz(channel);

  % Construct gammatone filter with one ERBaud bandwidth:
  gammafilter = ...
      Gfb_Filter_new(sampling_frequency_hz, center_frequency_hz, ...
                     gamma_order, bandwidth_factor);
  gammafilter.normalization_factor = ...
      gammafilter.normalization_factor / normalization_divisor;
  analyzer = Gfb_Analyzer_set_filter(analyzer, channel, gammafilter);
end


%%-----------------------------------------------------------------------------
%%
%%   Copyright (C) 2002, 2003  AG Medizinische Physik,
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
