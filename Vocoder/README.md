# Braecker-Vocoder #

### What is the Braecker-Vocoder? ###

The Braecker-Vocoder simulates different CI-coding strategies
and is able to either output an electrodogram or an auralized version 
of the electrodogramm, which normal hearing people can listen to via headphones.

### Supported coding strategies ###

Currently the following coding strategies are supported:

* CIS-coding strategy
* electro-acoustic stimulation (with CIS coding in the higher frequency bands)
* FSP-coding
* FAST like coding (maxima of envelope are sampled)
* FAST like coding with fine-structure sampling in the lower frequency bands

All coding strategies additionally can be used in an ACE like N-of-M electrode
selection process. For a complete list of all available coding strategies, look
at the function coding_strategy() in createCI.m.

Please note that all coding strategies are proof-of-concept coding strategies.
There is no gurantee, that they work exactly as the coding strategies available
in CI research interfaces or in CI speech processors. 

### Adjustable CI parameters ###

The most common CI parameter can be set:
 * Analysis frequency channels
 * M and T-levels
 * Number of electrodes
 * Number of active electrodes (n-of-m)
 * Total stimulation rate (vocoder sampling frequency)
 * pps per channel (only for CIS and FSP-coding strategy)
 * pulse-length (in samples)
 * Electrode-to-place mapping (auralistion only)
 * Simulation of spatial spread (auralisation only)

M and T-Levels can be defined per electrode. This will compress the signal
into CU-values. Note that there is currently no AGC, e.g. blockwise gain to
fit louder/softer sounds into the optimal audible range.
If you want to auralize a compressed sound, the CU-electrodogramm gets
inverted to normal amplitude. Because this is a nonlinear function, signal distortions
will occur and be audible. Of course, you can turn off the CU-conversion completely.
Note that you cannot use that electrodogramm for direct CI Stimulation.

To see all adjustable CI parameters, have a look at the function
setParameter() in createCI.m.

### How do I get started? ###

* Download the zip-File or clone the git repo locally.
* Open Matlab and cd to the downloaded directory

The most important file is `createCI.m`.
This function returns several function handles, which can be used to create
the CI simulation you want.  
Below are some code examples to get started with using the function:

Example Code for a complete Vocoder (CI simulation and auralisation)

```
#!matlab
[signal,fs] = wavread('OLSA.wav');
CI = createCI(); %Returns Function handels
vocoded_signal = CI.Vocoder(signal(:,1),fs,'fast_locked', 'debug',1);
sound(vocoded_signal,fs);
```

Example Code for producing only the electrodogramm:

```
#!matlab
[signal,fs] = wavread('OLSA.wav');
CI = createCI(); %Returns Function handels
parameter = CI.setParameter(signal(:,1),fs,'fast_locked','debug',1);
% in parameter are now all necessary infos for the CI simulation
[electrodogramm,parameter] = CI.Simulation(signal(:,1),fs,'fast_locked',parameter);
```
Now, we have an electrodogramm and we could auralize it:

```
#!matlab
vocoded_signal = CI.Auralisation(electrodogramm,parameter);
sound(vocoded_signal,parameter.voc_sampling_frequency_hz);
```

Example for creating binaural shifted electrodogramms:

```
#!matlab
[signal,fs] = wavread('OLSA.wav');
CI = createCI(); %Returns Function handels
[electrodogramm_l, electrodogramm_r] = CI.localizationEnhancement(signal(:,1),fs,'fast_locked', 45, 'debug', 1);
```

For further documentation, have a look at the [wiki](https://bitbucket.org/bwilliges/itd-sensitive-coding-strategy/wiki/Home).
### Contribution guidelines ###

* Whenever you change code, make sure the Testscript `Test_createCI.m` still runs.
* If you add an function, please also add a test-function in the Testscript. 
*  [Markdown Link](https://bitbucket.org/tutorials/markdowndemo)
* To be discussed