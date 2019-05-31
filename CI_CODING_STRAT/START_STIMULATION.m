%% Script to start CI coding strategies CIS and n-of-m
% Coding for MED EL 12 Electrode Implants
% Masterthesis Anja Eichenauer, University Oldenburg
% Date: 05.02.2016

%%
[Input, fs]=wavread('dahdlow');

%% Parameter allocation
Strategy=1; % 1 = CIS all electrodes on, 2= CIS electrodes 7 and 8 off, 3 = n-of-m 7-of-12, 4 = FL
T_Thresh=100*ones(22,1); % Threshold levels per Electrode [CU]
C_Thresh=1000*ones(22,1); % Comfortable levels per electrode [CU]
CalibLevel=10; % Base level + x dB [dB] (dynamic = 40 dB)
switch Strategy
    case 1
        % CIS all electrodes on, 900 pps
        par=CIParameters_1;
        par.C=C_Thresh;
        par.T=T_Thresh;
        par.vol=0.6;
    case 2
        % CIS electrodes 7 and 8 off, 500 pps
        par=CIParameters_2;
        par.C=C_Thresh;
        par.T=T_Thresh;
        par.vol=0.6;
    case 3
        % n-of-m 900 pps, 7-of-12
        par=CIParameters_3;
        par.C=C_Thresh;
        par.T=T_Thresh;
        par.vol=1;
    case 4
        % FL Strategy
        par=CIParameters_4;
        par.C=1000*ones(36,1);
        par.T=100*ones(36,1);
        par.vol=1;
end

%% Calibration: Set RMS-Level
DesiredLevel=20*log10(par.B)+CalibLevel;
desiredAmplitude=10^(DesiredLevel/20);
InputRMS=sqrt(mean(Input.^2)); % Root mean square
Input=Input/InputRMS*desiredAmplitude;

%% CI processing
[electrodogram, vDistance, StimOrder,vTime, mEnv]=CIStrat(Input,par,fs);
PlotElectrodogram(electrodogram, vTime)
%% Wrapper
% Send stimulation to RIB