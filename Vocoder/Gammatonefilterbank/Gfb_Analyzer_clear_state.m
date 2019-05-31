function analyzer = Gfb_Analyzer_clear_state(analyzer)
% analyzer = Gfb_Analyzer_clear_state(analyzer)
% Input (and then output) may also be a structure arrays of several pluggable
% analyzers.
%
% method for resetting the filter states to zeros
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002

% filename : Gfb_Analyzer_clear_state.m

for analyzer_index = 1:length(analyzer)
  for channel = [1:length(analyzer(analyzer_index).center_frequencies_hz)]
    filter = Gfb_Analyzer_get_filter(analyzer(analyzer_index), channel);
    filter = Gfb_Filter_clear_state(filter);
    analyzer(analyzer_index) = ...
      Gfb_Analyzer_set_filter(analyzer(analyzer_index), channel, filter);
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
