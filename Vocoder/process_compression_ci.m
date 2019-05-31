function C = process_compression_ci(A,B,M,alpha_c)

% usage:  C = process_compression_ci(A)
% input:  A = the amplitude values in each channel to simulate
%         B = Base level
%         M = saturation level
%         alpha = controls the steepness of the function
% output: C = clinical current unit
% 
% This function compresses the envelopes into the electrical dynamic range
% of CI users. 
%
% Reference: Laneau, 2005 PhD-thesis
%
% Author: Stefan Fredelake
% Date: 23-10-2008
% 

%THIS IS STEFANS CODE, INDEPENDENT OF FREQ-SPECIFIC TCL, MCL
% C = zeros(size(A));
% C(A<0.0156) = 0;                                           % eq. 5.17
% index = find(A>=0.0156 & A<0.5859);                              % eq. 5.17
% C(index) = log( 1 + 416.2 * ((A(index)-0.0156)/0.5703) )/6.033;  % eq. 5.17
% C(A>=0.5859) = 1;                                          % eq. 5.17

%Code from the dissertation of Brett Swanson (2008), p. 67 eq. 4.25
C = zeros(size(A));
C(A<B) = 0;
index = find(A>=B & A<M); 
C(index) = log( 1 + alpha_c * ((A(index)-B)/(M-B) ))/log(1+alpha_c);
C(A>=M) = 1; 

% Copyright (C) 2008 Stefan Fredelake, Oldenburg University
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 3 of the License, or (at your 
% option) any later version.
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
% You should have received a copy of the GNU General Public License 
% along with this program; if not, see http://www.gnu.org/licenses/.
% 