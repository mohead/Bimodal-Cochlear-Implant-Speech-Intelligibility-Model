function y = rms2(x, dim) 

% Usage: y = rms(x, dim) 
%
% for vector or matrix x, return its root mean square y,
% formula:  y = sqrt(mean(x.*conj(x)));
%
% second argument dim (necessary for 2-dim matrix input) 
% specifies the dimension along which to calculate the RMS
%
% Author:  Martin Hansen
% Date:    03 Nov 1999
% Updated: < 4 Apr 2005>
% Updated: < 6 Dec 2007>

%% This file is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published
%% by the Free Software Foundation; either version 2 of the License, 
%% or (at your option) any later version.
%% See the GNU General Public License for more details:
%% http://www.gnu.org/licenses/gpl

if min(size(x)) > 1,
  y = sqrt(mean(x.*conj(x), dim));  % works for complex vectors, too
else
  y = sqrt(mean(x.*conj(x)));       % works for complex vectors, too
end

if nargout < 1,
  for k=1:length(y),
    fprintf('*** info:   RMS = %f  (%g dB re. FS)\n', y(k), 20*log10(y(k)));
  end
end

% End of file:  rms.m
