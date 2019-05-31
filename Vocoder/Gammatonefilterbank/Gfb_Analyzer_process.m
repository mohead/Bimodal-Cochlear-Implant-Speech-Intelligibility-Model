function [output, analyzer] = Gfb_Analyzer_process(analyzer, input)
% [output, analyzer] = Gfb_Analyzer_process(analyzer, input)
%
% The analyzer processes the input data.
%
% PARAMETERS
% analyzer A Gfb_Analyzer struct created from Gfb_Analyzer_new. The
%          analyzer will be returned (with updated filter states) as
%          the second return parameter
% input   Either a row vector containing the input signal to process, or
%         a matrix containing different input signals for the different
%         channels.  Different rows correspond to different filter bands,
%         while different colums correspond to different instants of time. 
% output  A matrix containing the analyzer's output signals.  The
%         rows correspond to the filter channels
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002, Sep 2003

% filename : Gfb_Analyzer_process.m


if (analyzer.fast)
  % use a matlab/octave extension for fast computation.
  [output, analyzer] = Gfb_Analyzer_fprocess(analyzer, input);
else
  number_of_channels = length(analyzer.center_frequencies_hz);
  if (size(input,1) == 1)
    input_indizes = ones(1,number_of_channels);
  elseif (size(input,1) == number_of_channels)
    input_indizes = [1:number_of_channels];
  elseif (size(input,2) == 1)
    warning('Gfb:Analyzer:process', ...
            ['Input Signal is a column vector with ', ...
             'length != number_of_channels. Will convert it to row ' ...
             'vector']);
    input = input';
    input_indizes = ones(1,number_of_channels);
  else
    error('Gfb:Analyzer:process', ...
          ['Input Signal is a matrix with number of rows ~= number of ' ...
           'channels of this analyzer']);
  end
  output = zeros(number_of_channels, length(input));
  for channel = [1:number_of_channels]
    filter = Gfb_Analyzer_get_filter(analyzer, channel);
    [output(channel,:), filter] = ...
        Gfb_Filter_process(filter, input(input_indizes(channel),:));
    analyzer = Gfb_Analyzer_set_filter(analyzer, channel, filter);
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
