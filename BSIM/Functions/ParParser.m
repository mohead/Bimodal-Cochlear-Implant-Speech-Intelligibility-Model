function stControl = ParParser(sig_len,varargin)
% input parser for runBSIM:
p = inputParser;
%               Parameter name  , Default Value , Test method
p.addRequired('sig_len',@isnumeric); %signal length
p.addParameter('frame'          ,[sig_len sig_len 0]            ,@isnumeric);	%
p.addParameter('audiogram'      ,[0 -Inf -Inf; 20000 -Inf -Inf]	,@(x) (isnumeric(x)&&size(x,2)==3));	% audiogram in the form [frequency HL_left HL_right]
p.addParameter('CI_audiogram'   ,[0 -Inf -Inf; 20000 -Inf -Inf]	,@(x) (isnumeric(x)&&size(x,2)==3));	% audiogram in the form [frequency HL_left HL_right]
p.addParameter('display'        ,'text'         , @isstr);      % 
p.addParameter('errorflag'      ,true              , @islogical);	%
p.addParameter('ECprocessing'   ,true              , @islogical);	%
p.addParameter('BetterEar'      ,1              , @(x) ismember(x,[0 1 2 3]));	% 0 = no 1 = switch between left/right ear, 2 = left ear only, 3 right ear only
p.addParameter('angle'          ,0              , @isnumeric); 	%  
p.addParameter('short_term'     ,false              , @islogical);	%
p.addParameter('SIIparam'       ,[0 0.25]'       ,@(x) (size(x,2) == 2));	% Sii parameter value for left and right ear
p.addParameter('control'        ,[]             , @isstruct);	%
p.addParameter('subject_group'  ,'NH_bin'       , @isstr);      %
p.addParameter('plotfigures'    ,'none'         , @isstr);     	%
p.addParameter('plotfreqband'   ,8              , @(x) (isnumeric(x) && min(x) >= 0 && max(x) <= 30)); 	%
p.addParameter('plotangle'      ,0              , @isnumeric);	%
p.addParameter('impfkt'         ,2              , @isnumeric);	% Selected Bandimportance function of the SII
p.addParameter('window'         ,'hannwin'      , @isstr);      %
p.addParameter('bSkipBSIM'      ,false          , @islogical);  %
p.addParameter('bNewLevel'      ,false          , @islogical);  %
p.addParameter('bFrameParSet'   ,false          , @islogical);  %
p.addParameter('level'          ,[]             , @isscalar);   %
p.addParameter('fs'             ,44100          , @isscalar);	%
p.addParameter('framesec'       ,[]             , @isnumeric);	%
p.addParameter('data'           ,[]             , @isstruct);   %
p.addParameter('Bimodal_SII_Switch_value', 21, @isnumeric);
p.addParameter('gtf_FilterType' ,'gtf'          , @isstr);      %
p.addParameter('gtf_Gfb_cf'     ,[146.02 188.74 236.34 289.36 348.42 414.21 487.50 569.15 660.10 761.41 874.27 1000.00 1140.06 1296.07 1469.87 1663.48 1879.16 2119.41 2387.05 2685.19 3017.31 3387.29 3799.43 4258.55 4769.99 5339.72 5974.39 6681.39 7468.98 8346.32].'...
                                                , @isnumeric);  % default Gfb_center_frequencies(1/gtf_fERBBwRatio,140,1000,9000)
p.addParameter('gtf_fERBBwRatio',1              , @isstr);      %
p.parse(sig_len,varargin{:})
par = p; 
stControl                       = par.Results.control;
stControl.run.audiogram_QCKBSIM = par.Results.audiogram;
stControl.signal.iFs            = par.Results.fs;
stControl.model.EC.bErrorFlag   = par.Results.errorflag;
stControl.model.EC.EC_flank     = par.Results.ECprocessing;
stControl.model.EC.BetterEar    = par.Results.BetterEar;
stControl.model.sFilterType     = par.Results.gtf_FilterType;
stControl.model.FB.vfCenterFreqs= par.Results.gtf_Gfb_cf;
stControl.model.FB.fERBBwRatio  = par.Results.gtf_fERBBwRatio;
stControl.run.sDisplay          = par.Results.display;
stControl.run.angle             = par.Results.angle;          % angle is used for plots and BIF (nImpFkt)
stControl.run.subject_group     = par.Results.subject_group;
stControl.run.mfCI_audiogram    = par.Results.CI_audiogram;
stControl.run.freqBand          = par.Results.plotfreqband;
stControl.run.plotFigures       = par.Results.plotfigures;
stControl.run.plotAngle         = par.Results.plotangle;
stControl.run.short_term        = par.Results.short_term;
stControl.run.SIIval            = par.Results.SIIparam;
stControl.signal.iFrameLen                       = par.Results.frame(1);
stControl.signal.iFrameShift                     = par.Results.frame(2);

stControl.Bimodal_SII_Switch_value        = par.Results.Bimodal_SII_Switch_value;
stControl.model.ImpFkt          = par.Results.impfkt;
sWindowType                    	= par.Results.window;
stControl.run.bSkipBSIM         = par.Results.bSkipBSIM;  %
bNewLevel                       = par.Results.bNewLevel;  %

if ~isempty(par.Results.level)
    vfLevel      = par.Results.level(:);
    bNewLevel    = true;
end

idx = find(strcmp(varargin,'control'));
if ~isempty(idx) && idx+1 <= nargin && isstruct(varargin{idx+1})
    stControl = varargin{idx+1};
end

if ~isempty(par.Results.data)
    stControl    = par.Results.data.stControl;
    stResult     = par.Results.data.stResult;
    iNumFrames   = par.Results.data.iNumFrames;
    if ~bNewLevel
        vfLevel  = par.Results.data.vfLevel;
    end
    bSkipBSIM = true;
    iFrameLen = stControl.signal.iSigLen;
    iFrameShift = iFrameLen;
end 

if ~isempty(par.Results.framesec)
    if length(par.Results.frame(3))==1
        error('frame parameters already set! (double parameter?)');
    else 
        stControl.signal.iFrameLen   = round(par.Results.framesec(1) * stControl.signal.iFs);
        stControl.signal.iFrameShift = round(par.Results.framesec(2) * stControl.signal.iFs);
    end 
end 
stControl.HL.mfAudiogram = [0 -Inf -Inf; 20000 -Inf -Inf]; % audiogram used in BSIM03

% -------------------------------------------------------------------------
% constants dependent on input parameters
% -------------------------------------------------------------------------

if strcmp(sWindowType,'rectwin')
    vWindow = eval(sprintf('window(@%s,iFrameLen)',sWindowType));
    iFrameshift = iFrameLen;         % is this an error? iFrameShift instead?
else
    vWindow = hann(stControl.signal.iFrameLen,'periodic');
end
stControl.signal.vWindow = vWindow ./ sqrt(mean(abs(vWindow).^2));
end

