function y = hanwin(sig,flank_n) 

% window signal with hanning window - the complete window or with flanks
% ----------------------------------------------------------------------
% Usage: y = hanwin(sig)
%        y = hanwin(sig, flank_n)
%
%   input:  
%       sig      signal vector, or matrix, typically of size (N,1) or (N,2)
%       flank_n  (optional) length of hanning flank in samples.  
%                If flank_n > N/2, or if no flank_n is given, a
%                window with flank_n = N/2 is generated, i.e. a full
%                hanning window without "ones" in the middle is applied.
%  output:  
%       y        hanning-flanked windowed signal vector, or
%                matrix of hanning-flanked COLUMNS of sig 
%
% Copyright (C) 2003 by Martin Hansen, Fachhochschule OOW
% Author :  Martin Hansen <psylab AT jade-hs.de>
% Date   :  07 Oct 2003
% Updated:  <10 Apr 2005 19:38, mh>

%% This file is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published
%% by the Free Software Foundation; either version 2 of the License, 
%% or (at your option) any later version.
%% See the GNU General Public License for more details:
%% http://www.gnu.org/licenses/gpl

slen = size(sig);
% $$$ if slen(1) < slen(2) | slen >2,
% $$$   fprintf('error, hanwin geht nur für column! signal-vektoren, mit siglen(dim1) < siglen(dim2)\n');
% $$$   return;
% $$$ end
siglen=slen(1);
n_chan=slen(2);

if nargin < 2,
  flank_n = floor(siglen/2);
end

if 2*flank_n > siglen,
  fprintf(' *** warning, 2*flank_n > siglen, flanking now WHOLE signal\n');
  winlen = siglen;
  hwin = hanning_window(winlen);
else
  winlen = 2*flank_n;
  hflanks = hanning_window(winlen);
  hwin = [hflanks(1:flank_n); ones(siglen-winlen,1); hflanks(flank_n+1:end)];
end

% matricize
if n_chan>1,
  hwin = hwin*ones(1,n_chan);
end

y = sig.*hwin;  
  

function h = hanning_window(N)
% return N point symmetric hanning window including starting and
% ending zero

h = 0.5*(1+ cos(-pi+2*pi/(N-1)*(0:N-1)'));


% End of file:  hanwin.m

% Local Variables:
% time-stamp-pattern: "30/Updated:  <%2d %3b %:y %02H:%02M, %u>"
% End:
