function varargout = ...
    Gfb_Analyzer_split(analyzer, varargin)
% [analyzer1, analyzer2, ...] = ...
%     Gfb_Analyzer_split(analyzer, filter_order1 [, filter_order2, ... ] )
%
% Split a Gfb_Analyzer structure into several pluggable parts, so that
%
%    [outputN, analyzer] = Gfb_Analyzer_process(analyzer, signal);
% has the same effect as
%
%    [output1, analyzer1] = Gfb_Analyzer_process(analyzer1, signal);
%    [output2, analyzer2] = Gfb_Analyzer_process(analyzer2, output1);
%    ...
%    [outputN, analyzerN] = Gfb_Analyzer_process(analyzerN, output{N-1});
%
% while the intermediate signals also have the correct normalization.
% The filter states of the new anlyzers are zeroed out, however.
%
% PARAMETERS:
%     analyzer                   
%         The original Gfb_Analyzer object.
%     filter_orderX
%         The gammatone filter order for analyzerX.
%         If the sum of all given filter orders is less than the order of
%         the original analyzer, then an extra filter_order argument, taking up
%         the remaining orders, will be generated automatically.
% OUTPUT:
%     analyzer1
%         Gfb_Analyzer structure equivalent to the first filter_order1 orders
%         of all gammatone filters of analyzer.
%     analyzerX
%         Analyzers corresponding to the remaining filter orders. The analyzers
%         returned by this function are pluggable to each other if their order
%         is preserved.
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Sep Nov 2003

% filename : Gfb_Analyzer_split.m


orders = cell2mat(varargin);

filter = Gfb_Analyzer_get_filter(analyzer, 1);
%returns the sum of oders to the filter.gamma_order (e.g. 6)
filter.gamma_order = sum(orders);
if sum(orders) < filter.gamma_order
    orders = [orders, filter.gamma_order-sum(orders)];
elseif sum(orders) > filter.gamma_order
    error('Sum of given filter orders greater than order of original filterbank')
end

for index = 1:length(orders)
    varargout{index} = ...
        Gfb_Analyzer_reduce_to(analyzer, orders(index), index == 1);
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
