function Gfb_plot(figure_number, axis_vector, title_string, ...
		  xaxis_label, yaxis_label, domain, value)
% function Gfb_plot(figure_number, axis_vector, title_string, ...
%                   xaxis_labal, yaxis_label, domain, value)
%
% this function portably creates plots with labeled axes for matlab and octave
%
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan 2002

% filename : Gfb_plot.m

if (exist('OCTAVE_VERSION'))
    if (gnuplot_has_frames == 0)
        disp('WARNING: variable gnuplot_has_frames is 0. Assuming this is a');
        disp('         configuration error and setting it to 1.');
        gnuplot_has_frames = 1;
    end
    figure(figure_number);
    axis(axis_vector);
    title(title_string);
    xlabel(xaxis_label);
    ylabel(yaxis_label);
    plot(domain, value);
else
    figure(figure_number);
    plot(domain, value);
    axis(axis_vector);
    title(title_string);
    xlabel(xaxis_label);
    ylabel(yaxis_label);
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
