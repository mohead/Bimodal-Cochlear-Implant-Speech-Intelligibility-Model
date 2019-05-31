# Bimodal speech intelligibility model. 
This is a functional MATLAB model which estimates the speech reception threshold (SRT) of bimodal cochlear implant users [1]. It is based on the binaural speech intelligibility model (BSIM) [2]. This model predicts the SRT required in dB SNR required to achieve a certain percentage score on the Oldenburg sentence test (OlSa). It tries to achieve a reference speech intelligibility index (SII) value that is preset. The default values are set for 20%, 50%, and 80% correct answers on the OlSa. The model will use a sound signal that has the spectrum of speech as the noise and speech (Olnoise). It is also possible to use speech sentences as the speech signal, however, this would require averaging the results over several sentences. 

The model can predict the SRTs of Normal hearing (NH), aided hearing impaired (aHI), and vocoded cochlear implant for ACE (used in Cochlea models) and CIS (used in MED-EL) coding strategies on either ears. The model includes everything needed to replicate the results in [1]. 

## Getting started: 
Running createPlotData.m will recreate the results presented in [1]. To simulate a desired subject group call BimodalSIM(callname); Where call name has two parts, the mode of hearing and the specific case  separated by underscore symbols: M_M_S_S. The mode (M) of hearing can be one of the following: NH, HL (for aided hearing loss), CI and Mon (for monaural, i.e. deaf ear). The second part (S) can be NH, MHL (Mild Hearing loss), SHL (Severe hearing loss), Cochlear or MED-EL. The following are the examples simulated in [1]. 

        callname          = 'HL_Mon_SHL_NH';
Monaural aided severe hearing loss with the hearing loss on the left side.

        callname          = 'HL_Mon_MHL_NH';
Monaural aided mild hearing loss with the hearing loss on the left side.

        callname          = 'NH_Mon_NH_NH';
Monaural aided normal hearing with NH on the left side.

        callname          = 'Mon_NH_NH_NH';
Monaural aided normal hearing with NH on the right side.

        callname          = 'Mon_CI_NH_Cochlear';
Monaural CI (ACE) on the right side side.

        callname          = 'Mon_CI_NH_MED-EL';
Monaural CI (CIS) on the right side side.

        callname          = 'NH_NH_NH_NH';
Binaural normal hearing

        callname          = 'HL_CI_SHL_MED-EL';
Bimodal hearing with severe hearing loss on the left side, and CI (CIS) on the right side.

        callname          = 'HL_CI_SHL_Cochlear';
Bimodal hearing with severe hearing loss on the left side, and CI (ACE) on the right side.

        callname          = 'HL_CI_MHL_MED-EL';
Bimodal hearing with mild hearing loss on the left side, and CI (CIS) on the right side.

        callname          = 'HL_CI_MHL_Cochlear';
Bimodal hearing with mild hearing loss on the left side, and CI (ACE) on the right side.

        callname          = 'NH_CI_NH_MED-EL';
Bimodal hearing with normal hearing on the left side, and CI (CIS) on the right side.

        callname          = 'NH_CI_NH_Cochlear'; 
Bimodal hearing with normal hearing on the left side, and CI (ACE) on the right side.

More parameters can be adjusted while calling the function: 

|'subject_group'                               |Type          |Default value           |description| 
|----------------------------------------------|--------------|------------------------|----------------|
|'room'                                         |char          |'anechoic'              |Room-type, currently only anechoic is supported|
|'Display'                                      |char          |'text'                  |Visualization type:'text','notext','fig'|
|'plotFigures'                                  |char          |'none'                  |'all': all plots, 'SNR': SNR plots alone, 'EC': EC step plots alone, or 'none' |
|'plotFreqBand'                                 |numeric       | 8                      |Select a band index to plot between 1 and 30, or 0 for all, only if plotFreqBand is set to true|
|'noise_azim'                                   |numeric       |[-90,0,90]              |Angle of the noise|
|'plotAngle'                                    |numeric       |-90                     |Select an angle to plot, only if plotFreqBand is set to true, has to be a value from noise_azim|
|'Nchanns'                                      |numeric       | 2                      |Number of channels generated by the HRTF function. Some OpenMHA algorithms need 4 input channels|
|'Bimodal_SII_Switch_value'                     |numeric       | 21                     |The EA weighting function midpoint, explained in [1]. |
|'short_time_BSIM_flag':                        |logical       | false                  |Flag to use short-time stBSIM2010 (true) or do Batch Processing BSIM2010 (false)|
|'error_flag'                                   |logical       | true                   |Flag to use processing errors in EC mechanism (true) or not (false)|
|'Use_Shadow_filtering'                         |logical       | true                   |Flag to use shadow filtering, in MHA preprocessing|
|'Use_HL_Simulations'                           |logical       | false                  |Flag to use hearing loss simulation (true) or not (false)|
|'plotSRT'                                      |logical       | false                  |Plot SRT results|


For example, to run the model with different angles of noise arrival (i.e. -45°, 0°, and 45°), call the model as follows: BimodalSIM(callname,'noise_azim',[-45,0,45]).


## Replicating paper results and simulating other scenarios: 
At the time of this model's release, the shadowfiltering module is still not included in openMHA. Therefore, precalculated wav files (N000_NH.wav, N000_HL2.wav, etc) are used to regenerate the exact results found in [1] which used a shadowfiltering technique to separately process the noise and speech signals with the hearing aid processing. The shadow filtering is only available in a commercial version of the MHA. To simulate aided hearing impaired subgroups with different noise and speech angles, the openMHA can be used. In order to do so, clone the contents of openMHA and compile it in the openMHA folder and set the 'Use_Shadow_filtering' flag to false. Steps to compile the openMHA can be found in the readme file. Note that in this case, the results will slightly differ from those of [1]. 

The openMHA release 4.5.1 (2017-07-17) was tested with this model and can be obtained for free from www.openMHA.org.
The closed MHA version used is MHA-bimodal-4.5.1-x86_64-linux-gcc-5_N and can be purchased from https://www.hoertech.de. 

## License and Credits

If you want to use this code in your own work, please cite this software, using the following DOI: https://doi.org/10.5281/zenodo.3236176.

## References

[1] A. Zedan, B. Williges, and T. Jürgens. "Modeling Speech Intelligibility of Simulated Bimodal and Single-Sided Deaf Cochlear Implant Users." Acta Acustica United with Acustica 104, no. 5 (September 1, 2018): 918-21. https://doi.org/10.3813/AAA.919256.
[2] R. Beutelmann, T. Brand, and B. Kollmeier. "Revision, Extension, and Evaluation of a Binaural Speech Intelligibility Model." The Journal of the Acoustical Society of America 127, no. 4 (2010): 2479-97. https://doi.org/10.1121/1.3295575.
