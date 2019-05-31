% Script to calculate center frequencies for CI-Electrodes using
% Greenwoods 1990 function to relate frequency to place on
% Basiliarmembrane.

%
% Author : Ben Williges <ben.williges@uni-oldenburg.de>
%%=======================================================================%%

% The place on the basiliarmembrane is determined by information on the
% design of Cochlear Nucleus hybrid L24 Electrode:
% Length of electrode array: 16 mm
% Number of electrodes: 22, spread over 14,5 mm active length
%
clear; close; clc;
% constants for greenwoods function:
A = 165.4;
k = 0.88;
a = 0.06;
% Place parameters: Greenwoods function starts from the apex of the
% cochlea (0 mm):
% I Assume a mean cochlear duct length of 32 mm (see Cochlear Duct Length Plots in Med-El-Electrode Guide)), first electrode
% (most apical electrode) is located
% at aproximately 32 mm - 15,5 mm (subtract 0,5 mm to allow room for the soft
% tip). Last electrode would be located at 32mm - 15,5 mm + 14,5 mm (basal
% end)
%

%
% Assume linear distribution of electrodes (real distribution is different:
% Med-El: Number of Electrodes: 19, Active length = 23,1 mm, Contact
% Spacing = 2,1 mm -> 18*2,1 mm != 23,1 mm (Flex 28)

% Cochlear Contour-Advance-Elektrodenträger: Number of Electrodes: 22, Active length = 15 mm, Contact
% Spacing = 0.6818mm, Basal diameter: 0.8mm, Apical diameter: 0.5mm
x_cochlear = (32-18.5:15.0/(22-1):32-18.5+15.0);


% Med-El: Flex 24 (EAS-Electrode), Active length = 20,9, Array length =
% Insertion depth = 24 mm, thus the first contact is assumed to be at 23,5
% mm electrode array depth
x_medel_short = (32-23.5:20.9/(12-1):32-23.5+20.9); %12 Electrodes
x_medel = (32-23.5:20.9/(19-1):32-23.5+20.9); % 19 Electrodes

center_frequency_cochlear = A.*(10.^(x_cochlear.*a)-k); 
center_frequency_medel = A.*(10.^(x_medel.*a)-k);
short_center_frequency_medel =  round(A.*(10.^(x_medel_short.*a)-k));
