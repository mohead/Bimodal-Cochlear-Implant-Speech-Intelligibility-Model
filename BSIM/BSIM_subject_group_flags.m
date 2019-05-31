function para = BSIM_subject_group_flags(subject_group,Use_HL_Simulations)
% This function only takes the subject group name as an input with the
% following format: 'Ty_Ty_Au_Au'. Where Ty is the type which is one of the
% following ['NH,Mon,HL,CI']. And Au is the audiogram with the following
% possible values ['NH,Cochlear,MED-EL,SHL,MHL'] which stand for Normal hearing,
% Cochlear, Med-El, Severe HL, or Mild hearing loss respectively. 
% it outputs a structure that contains the parameters that depend on the
% group's chracteristics. 
% The parameters are as follows: 
% 
% % EC_flag: 1 for EC-processing activated, 0 for EC-processing
% dectivated. It is usually set to 0 for monaural cases because no binaural
% EC processing should happen (also saves computation time)
% 
% % better_ear_flag: 1 for better Ear, 2 for Left Ear, 3 for Right Ear.
% It is usually be set to the hearing ear in monaural cases.
% This option will only be effective when EC_flag is set to 0. 
% 
% % Nchanns The number of microphone channels used in openMHA, see
% convolveHRIR.m.
% Maximum number of channels is 6 and the order is as follows: 
% [left-front,right-front,left-rear,right-rear,left-in-ear,right-in-ear]
% 
% % sii_vals         [Left Ear;Right Ear]; 
% SII values are specified for both ears. The currently used SII values are
% 1. SII_NH that is used for normal hearing and hearing aid subject groups
% 2. SII_CI that is used for cochlear implant subject groups
% 3. SII_Deaf that is used for the deaf ear in all monaural cases.
% Each column in the SII_XX array contians SII values for 20% 50% and 80%
% words correct respectively. 
%
% Speech intelligibility indecies are the reference values that
% SRTFromSII.m tries to achieve by varying the noise level. The noise level
% added or subtracted is the SRT output of the BSIM Model.
%
% % Audiogram: The audiogram that will be used inside QuickBSIM.m. The 
% variable Audiogram is a string that will be used to generate the 
% frequencies at which the hearing loss (HL) is specified and the HL in dB
% for each ear.  
% The frequency dependend HL defines frequency bands, which should 
% not contribute to speech intelligibility.After the monaural and binaural
% (EC-process) speech and noise level has been calculated, the hearing 
% thresholds gets added to the noise levels, thus decreasing the SNR in 
% this frequency bands. Currently, this is is needed especially for HG1_mon
% and the Vocoders, where we have frequency bands, which are not audible, 
% when you listen to them, however, the binaural EC-process can still
% operate at -200 dB and increase the SNR in there.

SII_NH    = [0.193; 0.260; 0.327];         % Normal hearing SII is chosen such that it fits the median  SRT for 50% S0N0 of monaural NH listeners
SII_CI    = SII_NH +0.1560;              % SII_CI is chosen such that it fits the median SRT for 50% S0N0 of simulated monaural Cochlear Vocoder 

    %  NH;HL;CI;MN;

    Parser_Cell = strsplit(subject_group,'_');

    if strcmp(Parser_Cell{1},'Mon')||strcmp(Parser_Cell{2},'Mon')
            para.EC_flag            = false;                                 % No EC processing is possible when one ear is has on signal
            para.better_ear_flag    = 4-find(strcmp(Parser_Cell,'Mon')); % 3 - the monoral ear index will set the better ear to the hearing ear
                                                                         % i.e. if the left is set to mon, 4-1 = 3: use the right ear.
                                                                         % i.e. if the right is set to mon, 4-2 = 2: use the left ear.                                                        
    else
        para.better_ear_flag    = 1;
        para.EC_flag    = true;
    end 
    
    
    % Set the SII reference values
    para.sii_vals           = [SII_NH SII_NH];

    if strcmp(Parser_Cell{1},'CI')
        para.sii_vals(:,1) = SII_CI;
    end
    if strcmp(Parser_Cell{2},'CI')
        para.sii_vals(:,2) = SII_CI;
    end
    
    % Set the Audiogram
    LeftHearingLoss  = Parser_Cell{3};
    RightHearingLoss = Parser_Cell{4};
    if Use_HL_Simulations % In case of simulatoins, we dont need added noise
        if strcmp(LeftHearingLoss,'MHL')
            LeftHearingLoss = 'NH';
        end 
        if strcmp(LeftHearingLoss,'SHL')
            LeftHearingLoss = 'NH';
        end 
        if strcmp(LeftHearingLoss,'MHL')
            RightHearingLoss = 'NH';
        end 
        if strcmp(LeftHearingLoss,'SHL')
            RightHearingLoss = 'NH';
        end 
    else

    end 
    
	para.Audiogram   = PTAudiogram([LeftHearingLoss,'_',RightHearingLoss]);
    % Left for testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(subject_group,'test')
        para.sii_vals           = [SII_NH SII_NH];
        para.Audiogram      	= 'NH';
        para.Nchanns            = 2;
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
