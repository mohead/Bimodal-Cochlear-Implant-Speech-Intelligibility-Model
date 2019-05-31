 function [mDB, mF] = isothr_oetting(varargin);
%
% function [mDB, mF] = isothr(mF, iVer);
%
% returns the Iso Thresholds 'mDB' for frequencies
% 'mF'. (on default 'mDB' will be returned at
% available frequencies). 
%
% mF  : frequency scale (optional).
% iVer: defines Iso Threshold version.
%       iVer==1: ISO 226 (1985-05-01) extended at higher frequencies
%                by data from Klaus Bethge Dissertation.
%       iVer==2: same as iVer==1 for compatibility resons.
%       iVer==3: same as iVer==1 for compatibility resons.
%       iVer==4: reflects iso226:2003(E) prelim, as of Jun 2008
%                as well as "DIN ISO 28961:2011-10 DRAFT" 
%                as well as "DIN ISO 226:2011-10 DRAFT".
%       iVer==5: same as iVer==4
%       (default iVer==1, but produces warning and pause if not set)
%       iVer==7: zeros
%       iVer==8: iso389-8:2003 Table C.1 equalisation of HDA200 and convert
%                dB SPL to dB HL
% mDB : Threshold in DB SPL for tones (ISO226) 
%
% author / date : jens-e. appell / 3.95
% changes: jens-e. appell and jan rennies / 2008-06-26 added iVer option
%                  and data for iVer==4 (iso226:2003(E))
%          Tarik Siebe 2016-06-16 added iVer option
%                  and data for iVer==7 (iso389-8:2003 Table C.1)
%           Ben Williges 2016-08-19 renamed function to isothr_oetting, in
%           order to avoid name conflict with isothr.m distributed with
%           MHA (same Values (except for < 20 Hz), but different function input).
if nargin < 2,
	beep
	warning('Param iVer not passed. Set 1 to reproduce data before 2008-06-26, Set 4 to get iso226:2003(E)');
	warning('Param iVer == 1 will be used. Press key to continue...');
%	pause
end;

varargin= nargdef(varargin, [], 1);
mF		= varargin{1};
iVer	= varargin{2};

if iVer>0 & iVer<4,
	%% values from 20 Hz to 12500 Hz are taken from ISO 226 (1985-05-01)
	%% values at 14000 Hz and 15000 Hz are taken from ISO-Threshold-table
	%% in Klaus Bethges thesis. 
	%% values at 0 and 20000 Hz are not taken from ISO Threshold contour !!
	vThr	=     [100  74.3 65.0  56.3   48.4 41.7 35.5  29.8 25.1 20.7  16.8 13.8 11.2 8.9   7.2 6.0 5.0  4.4 4.2 3.8  2.6 1.0 -1.2 -3.6 -3.9 -1.1 6.6 15.3 16.4 11.6    16.0 24.1    70.0]';
	vsF	=1000*[0.01 0.02 0.025 0.0315 0.04 0.05 0.063 0.08 0.1  0.125 0.16 0.2  0.25 0.315 0.4 0.5 0.63 0.8 1.0 1.25 1.6 2.0 2.5  3.15 4.0  5.0  6.3 8.0  10.  12.5    14.0 15.0    20.0]';
elseif iVer>3 & iVer<6,
	%% data from iso226:2003(E) which is still preliminary but already
	%% commonly used. The date was received from Harm Remmers in June 2008
	vThr	= [78.5   68.7  59.5  51.1    44  37.5  31.5  26.5  22.1  17.9  14.4  11.4   8.6  6.2  4.4   3.0     2.2  2.4  3.5       1.7 -1.3 -4.2     -6.0  -5.4 -1.5       6.0 12.6       13.9   12.3];
	vsF	    = [20     25    31.5    40    50    63    80   100   125   160   200   250   315  400  500   630     800 1000 1250      1600 2000 2500      3150 4000 5000      6300 8000      10000  12500];
	%% values from "DIN ISO 226:2011-10 DRAFT" are the same
	%%vThr= [ 78.5  68.7  59.5  51.1    44  37.5  31.5  26.5  22.1  17.9  14.4  11.4   8.6  6.2  4.4     3     2.2  2.4  3.5       1.7 -1.3 -4.2        -6 -5.4 -1.5         6 12.6       13.9              12.3];
	%%vsF	= [   20    25  31.5    40    50    63    80   100   125   160   200   250   315  400  500   630     800 1000 1250      1600 2000 2500      3150 4000 5000      6300 8000      10000             12500];
	%% values from "DIN ISO 28961:2011-10 DRAFT" are the same after round
   %%vThr= [   79    69    60    51    44    38    32    27    22    18    14    11     9    6    4     3   2   2    2    4    2    2   -1   -4   -6   -6   -5   -2    4    6   13    14   14    13    13    12];
   %%vsF = [   20    25  31,5    40    50    63    80   100   125   160   200   250   315  400  500   630 750 800 1000 1250 1500 1600 2000 2500 3000 3150 4000 5000 6000 6300 8000 9000 10000 11200 12000 12500];
elseif iVer == 7  %use no correction (= zeros) if the incoming sound is given already in dB HL
    vsF	    = [ 125   160   200   250   315  400  500   630  800 1000 1250      1600 2000 2500      3150 4000 5000      6300 8000      ];
    vThr	= zeros(size(vsF));
	
elseif iVer == 8
    vThr	= [5 4.5 4.5 4.5 5 5.5 2.5 2.5 3 3.5 2 5.5 5 6 7 13 14.5 11 8.5];
	vsF	    = [ 125   160   200   250   315  400  500   630  800 1000 1250      1600 2000 2500      3150 4000 5000      6300 8000      ];

else
	error('Param iVer must be in the range of [1;4]');
end;

if isempty(mF),
	mDB	= vThr;
	mF		= vsF;
else
	lsF	= min(vsF);
	hsF	= max(vsF);
	vlb						= mF(:)>=lsF;
	vhb						= mF(:)<=hsF;
	mDB						= mF;
	mDB(find(vlb & vhb))	= interp1(log10(vsF),vThr,log10(mF(find(vlb & vhb))),'spline');
	mDB(find(~vlb))		= vThr(1);
	mDB(find(~vhb))		= vThr(end);
end;

return

%%-------------------------------------------------------------------------
%%
%% Copyright (C) 2012 by Jens-E. Appell, Fraunhofer Gesellschaft e.V.
%%
%% Permission to use, copy, and distribute this software/file and its
%% documentation for any purpose without permission by the author
%% is strictly forbidden.
%%
%% Permission to modify the software is granted, but not the right to
%% distribute the modified code.
%%
%% This software is provided "as is" without expressed or implied warranty.
%%
%%
%% AUTHOR CONTACT
%%
%% Dr. rer. nat. Jens-E. Appell
%% 
%% Fraunhofer Institute for Digital Media Technology
%% Haus des Hörens
%% Marie-Curie-Strasse 2
%% 26129 Oldenburg
%% Germany
%% 
%% jens.appell@idmt.fraunhofer.de
%%
%%-------------------------------------------------------------------------
