function createPlotData()
% This function simulates SRTs for all possible measurements conditions from the Vocoder-study
% in Zedan, et al, 2018, not yet published using the Bimodal-BSIM framework.
% By default it will simulate all 12 different subject groups. The
% different subject groups corresponding to a number from 1 to 12 are
% described in the subfunction get_subj_group().
% For each subject group, the result is plotted and saved to
% /Model_plots/Results_OlSa_bimodal_simulations/SRT_*.mat, where * is
% replaced by the subject group name.
% This files are then later used to plot the comparison between measured
% Vocoder results and Bimodal-BSIM modeled Results (using the function
% PlotSRT).

noiseAzim       = [-90;0;+90];
TargetSRTs      = [0.50]; % in probability
numSimAngles    = length(noiseAzim);
numSimSRTs      = length(TargetSRTs);

subject_groups  = 1:12; 
numSimulationsrunspergroup  = 1;
assert(numSimulationsrunspergroup>=1,'At least one run per group is needed');
% 125  Hz_db_HL  250  Hz_db_HL  500  Hz_db_HL     750 Hz_db_HL
% 1000 Hz_db_HL  1500 Hz_db_HL	2000 Hz_db_HL
% 4000 Hz_db_HL  6000 Hz_db_HL	8000 Hz_db_HL
numPoints       = numSimSRTs*numSimAngles;
%  Preallocate: rows: 3 noise azimus* 3 SRT target values* number of
%  simulatoins* 12 subject groups. Columns: Testlist+subject_group+azimus+SRT+HA-Side
% model_data      = zeros(numSimSRTs*numSimAngles*numSimulations*NumGroup,16);
w = waitbar(0,'Generating SRTs');                       % create a new waitbar, w with 0% progress
model_data2   = zeros(numSimSRTs*numSimAngles*numSimulationsrunspergroup*length(subject_groups),17);
for n_subjGroup = subject_groups
    model_data      = zeros(numSimSRTs*numSimAngles*numSimulationsrunspergroup,17);
	for n = 1:numSimulationsrunspergroup
       [callname,subj_group_name ]= get_subj_group(n_subjGroup);

        SRT      = BimodalSIM(callname,'Bimodal_SII_Switch_value',21);
        SRT      = SRT(:,2); % Select the 50% target SII.
       %  Preallocate: rows: 3 noise azimus* 3 SRT target values* number of
        %  simulatoins* 12 subject groups. Columns: Testlist+subject_group+azimus+SRT+HA-Side
        row_index   = ((1+numPoints*(n-1)):(n*numPoints));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fill data in the following column order: 
        % ID,  GroupID,  NoiseAzimus,  SRTinSNRbB,  %correct,  HA_Side,  Audiogram (11 columns)
        model_data(row_index,:) = [row_index.'+(n_subjGroup-1)*numSimulationsrunspergroup*length(row_index),...
        repmat(n_subjGroup,numPoints,1),   repmat(noiseAzim,numSimSRTs,1) ,...
        reshape(SRT,[],1),   rectpulse(TargetSRTs,numSimAngles),   zeros(numPoints,1),   zeros(numPoints,11)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        model_data2 ((n_subjGroup-1)*numPoints+n:...
            (n_subjGroup-1)*numPoints+n+numPoints-1,:) = model_data(row_index,:);
    end
        w = waitbar(n_subjGroup/length(subject_groups),w,['iteration: ', num2str(n_subjGroup), ' out of ', num2str(length(subject_groups))]);
        save(['SRT_Results' filesep 'SRT_' subj_group_name],'model_data');

        
end
save(['completeResultSet'],'model_data2')

delete(w);


function [callname,subj_group_name ]= get_subj_group(num)
switch num
    case 1 
        subj_group_name = 'HG1_mon';
        callname      	= 'HL_Mon_SHL_NH';
    case 2 
        subj_group_name = 'HG2_mon';
        callname      	= 'HL_Mon_MHL_NH';
    case 3 
        subj_group_name = 'NH_mon';
        callname      	= 'NH_Mon_NH_NH';
    case 4
        subj_group_name = 'Cochlear_mon';
        callname      	= 'Mon_CI_NH_Cochlear';
    case 5 
        subj_group_name = 'MED-EL_mon';
        callname      	= 'Mon_CI_NH_MED-EL';
    case 6
        subj_group_name = 'NH_bin';
        callname      	= 'NH_NH_NH_NH';
    case 7 
        subj_group_name = 'HG1_MED-EL';
        callname      	= 'HL_CI_SHL_MED-EL';
    case 8 
        subj_group_name = 'HG1_Cochlear';
        callname      	= 'HL_CI_SHL_Cochlear';
    case 9 
        subj_group_name = 'HG2_MED-EL';
        callname      	= 'HL_CI_MHL_MED-EL';
    case 10 
        subj_group_name = 'HG2_Cochlear';
        callname      	= 'HL_CI_MHL_Cochlear';
    case 11
        subj_group_name = 'MED-EL_SSD';
        callname      	= 'NH_CI_NH_MED-EL';  
    case 12
        subj_group_name = 'Cochlear_SSD';
        callname      	= 'NH_CI_NH_Cochlear';  
end
