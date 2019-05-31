function analyzer = Gfb_Analyzer_reduce_to(analyzer, new_order, isfirst)
% analyzer = Gfb_Analyzer_reduce_to(analyzer, new_order, isfirst)
%
% Makes a copy of (new_order) filter orders from this gammatone filterbank
% analyzer. Set isfirst to nonzero if the copy is the first in a set of
% pluggable analyzers, and to zero if it only sees other cokmplex gammatone
% filterbank output.
%
% These two cases need different normalization factors
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Nov 2003

% filename : Gfb_Analyzer_reduce_to.m


divisor = 2;
factor = 1;

if isfirst
    factor = divisor;
end

for filter_index = 1:length(analyzer.filters)
    filter = Gfb_Analyzer_get_filter(analyzer, filter_index);
    filter.normalization_factor = ...
        (filter.normalization_factor / divisor) ^ ...
        (new_order / filter.gamma_order) * factor;
    filter.gamma_order = new_order;
    analyzer = Gfb_Analyzer_set_filter(analyzer, filter_index, ...
                                       Gfb_Filter_clear_state(filter));
end

%%-----------------------------------------------------------------------------
%%
%%   Copyright (C) 2003  AG Medizinische Physik,
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
%%   Author: Tobias Peters (tpeters@uni-oldenburg.de)
%%
%%-----------------------------------------------------------------------------
