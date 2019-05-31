%% main function ----------------------------------------------------------
function [stResult,stControl] = BSIM03(mfSignal,IntNoi,stControl)
%
%  function [stResult,stControl] = BSIM03(mfSignal,stControl)
%
%  copyright (c) 2008 Rainer Beutelmann, 
%  Medizinische Physik, Universit�t Oldenburg
%
%  The BSIM function calculates the effective SNR increase due to binaural listening
% (BSIM) according to Beutelmann and Brand (2006) and Beutelmann et al.
% (2010)
%
%  Beutelmann R, Brand T (2006). "Prediction of speech intelligibility in
%  spatial noise and reverberation for normal-hearing and hearing-impaired
%  listeners", J. Acoust. Soc. Am., 120(1), 331-342.
%
%  Beutelmann, R., Brand, T. and Kollmeier, B. (2010). "Revision, extension,
%  and evaluation of a binaural speech intelligibility model (BSIM)", 
%  J. Acoust. Soc. Am., 127(4), 2479-2497. 

%  Version 01 (11/2013): new commentation added 
%  by Christopher Hauth and Thomas Brand

%  Version 02 (04/2014): graphical display of internal EC parameters is now
%  avaliable. An Instant_BSIM, calculating the EC parameters based on mixed
%  signals (speech and noise) is incorporated
%  by Christopher Hauth and Thomas Brand

%  Department for Medical Physics
%  Carl von Ossietzky (CvO) University Oldenburg

% ^ 
% | this is external help (type 'help BSIM') to display
% =========================================================================
% | this is internal documentation for power users
% v


% mfSignal is a (iSigLen x 4)-matrix and organized: [SpeechLeft SpeechRight NoiseLeft NoiseRight]
% 
% CAUTION: all variables containing data in the *frequency domain* have to use fftshift
% or obey the convention [negative frequencies, 0, positive frequencies] 
% (*other than MATLAB standard*)
%%
    % set all internal constants via assignin (function at the end of this file)
    BSIM_Sub_Constants;
    
    %% Preparation
    % ---------------------------------------------------------------------
    % prepare the signals, calculate spectra
    % ---------------------------------------------------------------------
    % do some zeropadding depending on maximum delay
    stControl.signal.iSigLen = size(mfSignal,1);
    
    iZeroPadAmount   = round(stControl.model.EC.fMaxDelay * stControl.signal.iFs);

    % calculate number of FFT points (minimum = 4096 (4*blocksize) -> for
    % block processing in short-time BSIM (stBSIM))
    stControl.signal.nFFT = 2^nextpow2(size(mfSignal,1)+2*iZeroPadAmount);
    
    % use fft length of 1024 for graphics [not active]
    stControl.signal.nFFT_Image = 1024;
    
    % minimum fft length is 4096
    if stControl.signal.nFFT < 4096
        stControl.signal.nFFT = 4096;
    end
        
    % init some buffer variables for graphical display of internal
    % parameters and gammatone filter output (things need to be changed, NOT URGENT)
    
    if any(strcmp(stControl.run.sDisplay,'fig'))
        delay_buffer = zeros(30,1);
        alpha_buffer = zeros(30,1);
        SNR_mon_left = zeros(30,1);
        SNR_mon_right= zeros(30,1);
        SNR_binaural = zeros(30,1);
        period       = zeros(30,1);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % enable internal noise that is added to external noise before EC
    % parameters estimation and calculation (rmsdB should fit calibration)
   
    Audiogram = [125 0 0; 250 0 0; 500 0 0; 1000 0 0; 2000 0 0; 3000 0 0; 4000 0 0; 6000 0 0; 8000 0 0];
     %Internal noise corresponding to hearing threshold 
    InternalNoiseLev = 10.*log10(HLintNoise(stControl.run.mfCI_audiogram,stControl.model.FB.vfCenterFreqs)); %this is intensity so check for which log stControl.run.mfCI_audiogram-10 stControl.run.mfCI_audiogram
 
 % calculate rms value of the noise
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate spectra of speech (left and right ear) and noise (left and right ear):
    mfSpectra         = fftshift(fft(mfSignal,stControl.signal.nFFT),1)/stControl.signal.iFs;
    
    mfSpectraIntNoise = fftshift(fft(IntNoi,stControl.signal.nFFT),1)/stControl.signal.iFs;
    
    % internal parameter display
    if any(strcmp(stControl.run.sDisplay,'text'))
        disp(sprintf(' #:  Sp-Lv  Ns-Lv  mon-L  mon-R  SNR     alpha     tau     #search'));
    end
    
    stControl.band.bSaveNow = false;
 

%% main processing loop starts here:   
    % ----------------------------------------------------------
    % cycle through all frequency bands given in stControl
    % ----------------------------------------------------------
    for iBandCount = 1:length(stControl.model.FB.vfCenterFreqs)
        
        stControl.band.iCfIdx = iBandCount; % index of centre frequencies
        % calculate periodicity of corresponding centre frequency (used to find correct position of maximum SNR in EC-mechanism)
        stControl.band.EstPer = 1/stControl.model.FB.vfCenterFreqs(stControl.band.iCfIdx);
        period(iBandCount) = stControl.band.EstPer;
        % ----------------------------------------------------------
        % initialize
        % ----------------------------------------------------------
        bDone     = false;
        iCount    = 0;
        bPreCalc  = false;
        % ----------------------------------------------------------
        % build narrow band filter
        % ----------------------------------------------------------
        % get gammatone/third-octave band/critical band filter  
        [stControl, spvfTransferFct, spvfOmega] = BSIM_Sub_TransferFunction(stControl,stControl.model.FB.fERBBwRatio);
        viNZIdx = find(spvfOmega);
        
        % if graphical display, do the the same for a fft-length of 1024
        % samples
        
        if stControl.signal.nFFT < 4096
            stControl.signal.nFFT = 4096;
        end
        
        % -----------------------------------------------------------------
        % control centre frequency (plot transfer function of gammatone filter)
        % -----------------------------------------------------------------
%         figure(2) 
%         semilogx(spvfOmega/(2*pi),10*log10(abs(spvfTransferFct)));
%         hold on;
%         drawnow;
        
        % -----------------------------------------------------------------
        % calculate all overall parameters
        % -----------------------------------------------------------------
        % cut out the current frequency band from the full spectra
        % (apply gammatone filter)
        % matrix will be sparse
        % speech and noise separately
        
        spmfBandSpectra         = sparse(repmat(viNZIdx,1,4),repmat(1:4,size(viNZIdx,1),1),repmat(spvfTransferFct(viNZIdx),1,4) .* mfSpectra(viNZIdx,:),size(spvfTransferFct,1),4);
        spmfBandSpectraIntNoise = sparse(repmat(viNZIdx,1,2),repmat(1:2,size(viNZIdx,1),1),repmat(spvfTransferFct(viNZIdx),1,2) .* mfSpectraIntNoise(viNZIdx,:),size(spvfTransferFct,1),2);
        spmfBandSpectraIntNoise(:,1) = spmfBandSpectraIntNoise(:,1).*10.^((-65+InternalNoiseLev(iBandCount,1))/20);
        spmfBandSpectraIntNoise(:,2) = spmfBandSpectraIntNoise(:,2).*10.^((-65+InternalNoiseLev(iBandCount,2))/20);
        
        
        spmfBandSpectra(:,[3 4]) = spmfBandSpectra(:,[3 4])+ spmfBandSpectraIntNoise;
        % calculate overall intensities of all input  channels:
        vfIntensities = full(sum(abs(spmfBandSpectra).^2,1)*stControl.signal.iFs/stControl.signal.nFFT)+eps;  
        
        
        [stControl, vfInternalNoise] = BSIM_Sub_InternalNoise(stControl);        % get internal noise level from individual audiogram and internal noise spectrum level for normal hearing people
        vfIntensities([3 4]) = vfIntensities([3 4]) + vfInternalNoise;           % add internal noise level to external noise level  
                                                                                 % Operation is done for convenience here. In fact, the internal noise is set to zero (Hearing threshold to -Inf)
                                                                                 % for calculation of equalization parameters
        % Calculate Crosscorrelation of Speech [1 2](S_L(w).* conj(S_R(w))
        % see equation (11) in Beutelmann, R., Brand, T. and Kollmeier, B.
        % (2010), and the same for noise [3 4] 
        spmfXSpec = spmfBandSpectra(:,[1 3]) .* conj(spmfBandSpectra(:,[2 4]));
        
        
        % if there already exist values for alpha and tau from a previous calculation: 
        % run EC-mechanism with these parameters
        if isfield(stControl.model.EC,'Result') && isfield(stControl.model.EC.Result,'Tau') && isfield(stControl.model.EC.Result,'Alpha')
            
            fDelay = stControl.model.EC.Result.Tau(stControl.band.iCfIdx);
            tau    = fDelay;
            fAlpha = 10.^(stControl.model.EC.Result.Alpha(stControl.band.iCfIdx)/20);
            stControl.model.EC.bErrOn = true;
            % run EC-mechanism
            [stControl, speech_lev, noise_lev, tau, fAlpha] = BSIM_Sub_ECLevel(stControl,struct('bUseFFT',0,'bOptAlpha',0,'bMonoLimit',0),vfIntensities,fAlpha,tau,spvfOmega,spmfXSpec);
            stControl.model.EC.bErrOn = false;
            % compute SNR
            snr = speech_lev./noise_lev;
            c = 0;
            bPreCalc = true;

        % otherwise, estimate a value for tau (delay) and calculate the resulting
        % alpha:
        else
            % ----------------------------------------------------------
            % find a first coarse estimate of the snr maximum position
            % ----------------------------------------------------------
            [stControl, fAlpha, fDelay, fMaxSNR] = BSIM_Sub_CoarseEstMax(stControl,spmfXSpec,vfIntensities,spvfOmega);

            % ----------------------------------------------------------
            % search for the exact location of the maximum snr by quadratic
            % interpolation
            % ----------------------------------------------------------
          
            % indices 1 2 3 -> Delaypos(1 2 3)-> 2 = vfTau (needed for interpolation)
            viNewIdx = [1 2 3].';
            vfNewX = fDelay(:);
            vfNewY = 10*log10(fMaxSNR(:));

            iMaxCount = 3;
            c = 0;
            % ---------------------------------------------------

            while ~bDone

                if ~isempty(viNewIdx)
                    tau        = reshape(fDelay(viNewIdx(:)),1,[]);
                    iCount = iCount + 1;
                    
                    [stControl, speech_lev, noise_lev, tau, fAlpha] ...
                        = BSIM_Sub_ECLevel(stControl,struct('bUseFFT',0,'bOptAlpha',1,'bMonoLimit',0),...
                                           vfIntensities,[],tau,spvfOmega,spmfXSpec);
 
                    snr = speech_lev./noise_lev;
                    vfNewX(viNewIdx(:),1) = tau(:);
                    vfNewY(viNewIdx(:),1) = full(10*log10(snr(:)));
                    c = c + max(size(tau));
                end
                
                if any(isnan(vfNewX))
                    disp('vfNewX NaN');
                end
                if any(isnan(vfNewY))
                    disp('vfNewY NaN');
                end
% Intersample Interpolation by Quadratic approximation  (C. Implementation of BSIM, Beutelmann, R.,
% Brand, T. and Kollmeier, B.(2010))
                [stControl, vfNewX, vfNewY, viNewIdx, fCurve] = BSIM_Sub_QuadrOpt(stControl, vfNewX, vfNewY);

                bDone = ...
                    (fCurve < 0 && max(vfNewY(:))-min(vfNewY(:)) < 0.1) ...    % a maximum is found and the snrs differ only little
                    && ...
                    max(vfNewX)-min(vfNewX) < (pi/max(spvfOmega)) ...          % ... or the tau positions differ only little
                    || ...
                    iCount > iMaxCount ...                                     % ... or the maxmium number of iterations is reached
                    || ...
                    max(vfNewY-10*log10(max(vfIntensities(1:2)./vfIntensities(3:4)))) < 1;  % ... or none of the snrs is greater than the best monaural snr

                fDelay = vfNewX(:);
                vfNewY = vfNewY(:);

            end
            
            % Equalization-Cancellation Stage
            [stControl, speech_lev, noise_lev, vfNewX, fAlpha] = BSIM_Sub_ECLevel(stControl,struct('bUseFFT',0,'bOptAlpha',1,'bMonoLimit',0),vfIntensities,[],fDelay.',spvfOmega,spmfXSpec);
            if stControl.model.EC.bErrorFlag
            stControl.model.EC.bErrOn = true; %true
            end
            [stControl, speech_lev, noise_lev, vfNewX, fAlpha] = BSIM_Sub_ECLevel(stControl,struct('bUseFFT',0,'bOptAlpha',0,'bMonoLimit',0),vfIntensities,fAlpha,fDelay.',spvfOmega,spmfXSpec);
            stControl.model.EC.bErrOn = false;

        end
        snr_new = speech_lev/noise_lev;
        % result calculation 
        % ------------------
        fScaleFactor = stControl.signal.iFs/stControl.signal.iSigLen/stControl.model.FB.vfBandWidth(stControl.band.iCfIdx);
        fAlpha = full(fAlpha);
        snr   = full(snr);
        snr_new = full(snr_new);
        % monaural SNR
        mono  = [vfIntensities(1)/vfIntensities(3),vfIntensities(2)/vfIntensities(4)];
       
        [fMaxMonoSNR,iMaxMonoIdx] = max(mono);
        [stResult.SNR(stControl.band.iCfIdx,1),viBestSNRIdx] = max(snr);
        SNR_binaural(iBandCount) = 10*log10(max(max(snr_new)));
        
        
        % use monaural SNR if it is larger compared to the "binaural" SNR
        % (after EC-processing)
         if (fMaxMonoSNR < max(snr(:)))
             stResult.SpeechLev(stControl.band.iCfIdx,1) = (10*log10(speech_lev(viBestSNRIdx)*fScaleFactor));
             stResult.NoiseLev( stControl.band.iCfIdx,1) = 10*log10(noise_lev(viBestSNRIdx)*fScaleFactor);
         else
            stResult.SpeechLev(stControl.band.iCfIdx,1) = (10*log10(vfIntensities(iMaxMonoIdx  )*fScaleFactor));
            stResult.NoiseLev( stControl.band.iCfIdx,1) = (10*log10(vfIntensities(iMaxMonoIdx+2)*fScaleFactor));
         end

     
        stResult.SNR(stControl.band.iCfIdx,1) = stResult.SpeechLev(stControl.band.iCfIdx)-stResult.NoiseLev(stControl.band.iCfIdx);
        stResult.Mono( stControl.band.iCfIdx,:) = full(10*log10(mono));
        SNR_mon_left(iBandCount) = full(10*log10(mono(1)));
        SNR_mon_right(iBandCount)= full(10*log10(mono(2)));
        
        
        stResult.Alpha(stControl.band.iCfIdx,1) = dB(full(fAlpha(viBestSNRIdx)));
        stResult.Tau(stControl.band.iCfIdx,1) = full(fDelay(viBestSNRIdx));
        
       
        alpha_buffer(iBandCount) = dB(full(fAlpha(viBestSNRIdx)));
        delay_buffer(iBandCount) = full(fDelay(viBestSNRIdx));
       
        stControl.model.Delays(iBandCount,stControl.model.Counter) =  full(fDelay(viBestSNRIdx));
        stControl.model.Alphas(iBandCount,stControl.model.Counter) = dB(full(fAlpha(viBestSNRIdx)));        


        stControl.band.bSaveNow = true;
        
        if stControl.model.EC.bErrorFlag
        stControl.model.EC.bErrOn = true; %true
        end
        
        [stControl, sp_test, ns_test, tau_test, alpha_test] = ...
            BSIM_Sub_ECLevel(stControl,struct('bUseFFT',0,'bOptAlpha',0,'bMonoLimit',0),...
                             vfIntensities,fAlpha(viBestSNRIdx),stResult.Tau(stControl.band.iCfIdx),spvfOmega,spmfXSpec);
        stControl.model.EC.bErrOn = false;

        stResult.vfInt(stControl.band.iCfIdx,:)   = stControl.band.vfInt;
        stResult.vfCross(stControl.band.iCfIdx,:) = stControl.band.vfCross;
        
        
%% plot Cancelation step. 
    if ( (strcmp(stControl.run.plotFigures,'all')||strcmp(stControl.run.plotFigures,'EC') ) &&...
       (stControl.run.freqBand==iBandCount||stControl.run.freqBand==0) &&...
       (stControl.run.plotAngle==stControl.run.angle) ) 

        % Transfer the signals back to time domain
        mtSignal        = ifft(ifftshift(full(spmfBandSpectra),1),stControl.signal.nFFT).*stControl.signal.iFs;
        % Get the delay as an integer
        idelay          = -round(stControl.signal.iFs*fDelay(viBestSNRIdx));
        % Padd the signals in order to realize the delay
        mtRealSignal	= real(mtSignal);
        % Shift the left ear signal by the delay
        mtSigShifted    = [circshift(mtRealSignal(:,1),idelay), mtRealSignal(:,2), ...
                           circshift(mtRealSignal(:,3),idelay), mtRealSignal(:,4)];
        % Multiply the signal by the attentuation factor
        mtSigScaled     = [fAlpha(viBestSNRIdx).*mtSigShifted(:,1), mtSigShifted(:,2), ...
                           fAlpha(viBestSNRIdx).*mtSigShifted(:,3), mtSigShifted(:,4)];
        % Subtract the left ear from the right ear
        ecSig           = mtSigScaled(:,2)-mtSigScaled(:,1);
        ecNoi           = mtSigScaled(:,4)-mtSigScaled(:,3);
        % Plot the signal
        
        % Get maximum and minimum data and set data range. 

        % Get 2 periods of the signal for display
        dispSigLength   = round(2*stControl.band.EstPer*stControl.signal.iFs);
        % Calculate the x-axis limit values in order to show the difference
            % Calculate the starting limit
        XaxisMin        = round(stControl.signal.iSigLen/4)-round(dispSigLength/2);
            % Calculate the end limit
        XaxisMax        = XaxisMin+ dispSigLength;
        % Get 110% of the Y-axis range
        YaxisMax        = max(max(mtRealSignal(XaxisMin:XaxisMax,:)))*1.1;
        YaxisMin        = min(min(mtRealSignal(XaxisMin:XaxisMax,:)))*1.1;
        % Scale the x-axis limits to match the time of signal
        XaxisMin        = XaxisMin./stControl.signal.iFs;
        XaxisMax        = XaxisMax./stControl.signal.iFs;
        % Calculate the x-axis in time for plotting, and take the padding
        % length into account
        xAxis           = linspace(0,stControl.signal.nFFT/stControl.signal.iFs,stControl.signal.nFFT);

        figure('Name',['Subject Group: ' stControl.run.subject_group ', Angle: ' num2str(stControl.run.angle) ', Center Frequency ' num2str(stControl.model.FB.vfCenterFreqs(stControl.band.iCfIdx))]);

        % Create the subplots. 
        subplot(2,5,1)
        plot(xAxis,mtRealSignal(:,[1,2]))
        axis([XaxisMin XaxisMax YaxisMin YaxisMax])
        title('Real Part of Original Signal')
        legend('Left Ear Signal','Right Ear Signal')
        xlabel('Time')
        ylabel('Amplitude')
        grid on
        
        subplot(2,5,6)
        plot(xAxis,mtRealSignal(:,[3,4]))        
        axis([XaxisMin XaxisMax YaxisMin YaxisMax])
        title('Real Part of Original Noise')
        legend('Left Ear Noise','Right Ear Noise')
        xlabel('Time')
        ylabel('Amplitude')
        grid on
        
        subplot(2,5,2)
        plot(xAxis,mtSigShifted(:,[1,2]))        
        axis([XaxisMin XaxisMax YaxisMin YaxisMax])
        title('After Shift')
        legend('Left Ear Signal','Right Ear Signal')
        xlabel('Time')
        ylabel('Amplitude')
        grid on

        subplot(2,5,7)
        plot(xAxis,mtSigShifted(:,[3,4]))        
        axis([XaxisMin XaxisMax YaxisMin YaxisMax])
        title('After Shift')
        legend('Left Ear Noise','Right Ear Noise')
        xlabel('Time')
        ylabel('Amplitude')
        grid on
        
        subplot(2,5,3)
        plot(xAxis,mtSigScaled(:,[1,2]))       
        axis([XaxisMin XaxisMax YaxisMin YaxisMax])
        title('After Scaling')
        legend('Left Ear Signal','Right Ear Signal')
        xlabel('Time')
        ylabel('Amplitude')
        grid on
        
        subplot(2,5,8)
        plot(xAxis,mtSigScaled(:,[3,4]))    
        title('After Scaling')
        legend('Left Ear Noise','Right Ear Noise')
        axis([XaxisMin XaxisMax YaxisMin YaxisMax])
        xlabel('Time')
        ylabel('Amplitude')
        grid on
        
        subplot(2,5,4)
        plot(xAxis,mtSigScaled(:,2)-mtSigScaled(:,1))       
        axis([XaxisMin XaxisMax YaxisMin YaxisMax])
        title('Cancelation')
        legend('Left Ear Signal - Right Ear Noise')
        xlabel('Time')
        ylabel('Amplitude')
        grid on
        
        subplot(2,5,9)
        plot(xAxis,mtSigScaled(:,4)-mtSigScaled(:,3))       
        axis([XaxisMin XaxisMax YaxisMin YaxisMax])
        title('Cancelation')
        legend('Left Ear Signal - Right Ear Noise')
        xlabel('Time')
        ylabel('Amplitude')
        grid on
        
        subplot(2,5,5)
        plot(xAxis,mtRealSignal)
        title(['Original Signal at Central Freq ' num2str(stControl.model.FB.vfCenterFreqs(stControl.band.iCfIdx)) 'Hz'])
        legend('Left Ear Signal','Right Ear Signal','Left Ear Noise','Right Ear Noise')
        xlabel('Time')
        ylabel('Amplitude')
        grid on
        
        subplot(2,5,10)
        plot(xAxis,ecSig)
        hold on 
        plot(xAxis,ecNoi)
        hold off
        title('Equalized by Cancelation Signal')
        legend('Signal','Noise')
        xlabel('Time')
        ylabel('Amplitude')
        grid on


% fAlpha(viBestSNRIdx)
% fDelay(viBestSNRIdx)

end
%%  Plots        
        if any(strcmp(stControl.run.sDisplay,'text'))&&...
          (strcmp(stControl.run.plotFigures,'all')||strcmp(stControl.run.plotFigures,'SNR'))
            % print some useful information
            if (fMaxMonoSNR < max(snr(:)))
                disp(sprintf('%2u:  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f   %6.1f  %8.1f   %3u  %7.2f-%7.2f',...
                stControl.band.iCfIdx,...
                stResult.SpeechLev(stControl.band.iCfIdx  ),...
                stResult.NoiseLev( stControl.band.iCfIdx  ),...
                stResult.Mono(     stControl.band.iCfIdx,:),...
                stResult.SNR(      stControl.band.iCfIdx  ),...
                stResult.Alpha(    stControl.band.iCfIdx  ),...
                stResult.Tau(      stControl.band.iCfIdx  )/1e-6,...
                c,full([min(spvfOmega(viNZIdx)),max(spvfOmega(viNZIdx))]/2/pi)));
            else
                disp(sprintf('%2u:  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f  (%6.1f  %8.1f)  %3u  %7.2f-%7.2f',...
                stControl.band.iCfIdx,...
                stResult.SpeechLev(stControl.band.iCfIdx  ),...
                stResult.NoiseLev( stControl.band.iCfIdx  ),...
                stResult.Mono(     stControl.band.iCfIdx,:),...
                stResult.SNR(      stControl.band.iCfIdx  ),...
                stResult.Alpha(    stControl.band.iCfIdx  ),...
                stResult.Tau(      stControl.band.iCfIdx  )/1e-6,...
                c,full([min(spvfOmega(viNZIdx)),max(spvfOmega(viNZIdx))]/2/pi)));
            end
        end
    end  % end of main loop
    %
   % close(h);
    
    if any(strcmp(stControl.run.sDisplay,'fig'))&&...
    (strcmp(stControl.run.plotFigures,'all')||strcmp(stControl.run.plotFigures,'SNR'))

        angle = stControl.run.angle;
            short_term = stControl.run.short_term;
            
            if short_term==0
                H = figure();
            end
          
            
            if short_term==1
                
                 if angle<0;
                 H = figure(abs(angle));
                 else
                 H = figure(angle+1);
                 end
            end
            set(H,'name',[sprintf('angle=%d',angle) '  ' '''' stControl.run.subject_group '''']);
            
            % left ear channel:
            set(H,'units','normalized','Position',[0.01 -0.01 0.8 0.7]);
     
        % binaural part:
        bin1 = subplot(2,3,1);
        bin2 = subplot(2,3,2);
        bin3 = subplot(2,3,3);   
        bin4 = subplot(2,3,4);   
        
        set(bin1,'Position',[0.06 0.580 0.39 0.34]);        
        set(bin2,'Position',[0.06 0.100 0.39 0.34]);
        set(bin3,'Position',[0.55 0.580 0.39 0.34]);
        set(bin4,'Position',[0.55 0.100 0.39 0.34]);
       
        Delayplot = stControl.model.Delays(:,1:stControl.model.Counter)./(10^-6);
        
        if short_term ==1
 
          errorbar(stControl.model.FB.vfCenterFreqs,mean(Delayplot,2),std(Delayplot,1,2),'kx','Markersize',6,'Parent',bin1); 
          set(bin1,'NextPlot','add','YGrid','on','XGrid','on','XScale','log')
         % fill([stControl.model.FB.vfCenterFreqs;flipud(stControl.model.FB.vfCenterFreqs)],[mean(Delayplot,2)-2*std(Delayplot,1,2); flipud((mean(Delayplot,2) + 2*std(Delayplot,1,2)))], ...
         % [.5, .9 , 0],'FaceAlpha',0.1,'Parent',bin1);    
           set(bin1,'NextPlot','add','YGrid','on','XGrid','on','XScale','log')
        else
            set(bin1,'NextPlot','add','YGrid','on','XGrid','on','XScale','log');
            semilogx(stControl.model.FB.vfCenterFreqs,delay_buffer./(10^-6),'ko','Parent',bin1,'Markersize',6)
        end
     
        % plot period of gammatone filter as limit
        semilogx(bin1,stControl.model.FB.vfCenterFreqs(5:end),period(5:end)/2./(10^-6),'k','Linewidth',2);
        semilogx(bin1,stControl.model.FB.vfCenterFreqs(5:end),-period(5:end)/2./(10^-6),'k','Linewidth',2);
        ylabel(bin1,'\bf{Delay in \mu s}');
        title(bin1, sprintf('Delay EC stage; Noise source @ %d degrees',angle));%,mean(angle_est(1:15,:))))%f;
        ti=get(bin1,'title');
        set(ti,'FontWeight','bold')
        % alpha
        if short_term==1
            errorbar(stControl.model.FB.vfCenterFreqs,mean(stControl.model.Alphas(:,1:stControl.model.Counter),2),std(stControl.model.Alphas(:,1:stControl.model.Counter),1,2),'kx','Markersize',6,'Parent',bin2)
             set(bin2,'NextPlot','add','YGrid','on','XGrid','on','XScale','log')
         %  fill([stControl.model.FB.vfCenterFreqs;flipud(stControl.model.FB.vfCenterFreqs)],[mean(stControl.model.Alphas(:,1:stControl.model.Counter),2)-2*std(stControl.model.Alphas(:,1:stControl.model.Counter),1,2); flipud((mean(stControl.model.Alphas(:,1:stControl.model.Counter),2) + 2*std(stControl.model.Alphas(:,1:stControl.model.Counter),1,2)))], ...
          % [0.9, .4, 0],'FaceAlpha',0.1,'Parent',bin2);
        else
            set(bin2,'YGrid','on','XGrid','on','NextPlot','add','XScale','log');
            semilogx(bin2,stControl.model.FB.vfCenterFreqs,alpha_buffer,'ko','Markersize',6);
        end
        ylabel(bin2,'\bf{Attenuation in dB}');
        title (bin2,'\bf{Attenuation EC stage}');
        
        % SNRs
        semilogx(bin3,stControl.model.FB.vfCenterFreqs,SNR_binaural,'Color',[0 0 0],'Linewidth',2);
        set(bin3,'NextPlot','add','YGrid','on','XGrid','on')
        semilogx(bin3,stControl.model.FB.vfCenterFreqs,SNR_mon_left,'b','Linewidth',2);
        semilogx(bin3,stControl.model.FB.vfCenterFreqs,SNR_mon_right,'r','Linewidth',2);
        legend(bin3,'Bin','MonL','MonR','Location','Northoutside','Orientation','Horizontal');
        xlabel(bin3,'\bf{Frequency in Hz}');  
        ylabel(bin3,'\bf{SNR in dB}');
        title(bin3,'\bf{Signal-to-Noise Ratios}');
        
        
        % Binaural SNR Benefit
        SNR_bin_benefit = SNR_binaural-max(SNR_mon_left,SNR_mon_right);
        semilogx(bin4,stControl.model.FB.vfCenterFreqs,SNR_bin_benefit,'Color',[0 0 0],'Linewidth',2);
        set(bin4,'NextPlot','add','YGrid','on','XGrid','on')
        SNR_no_bin_ben  = SNR_bin_benefit;
        SNR_no_bin_ben(SNR_no_bin_ben>0)=nan;
        semilogx(bin4,stControl.model.FB.vfCenterFreqs,SNR_no_bin_ben,'r','Linewidth',2);
        legend(bin4,'Benefit','Hindrance','Location','Northoutside','Orientation','Horizontal');
        xlabel(bin4,'\bf{Frequency in Hz}');  
        ylabel(bin4,'\bf{SNR in dB}');
        title(bin4,'\bf{SNR Benefit of Binaural vs Monoral Hearing}');

        drawnow;
        stControl.model.Counter = stControl.model.Counter+1;
        %savefig(['..' filesep 'Figures' filesep stControl.run.subject_group '_' int2str(angle) '.fig']);
    end
    %%
% All Subroutines start here:

%% BSIM_Sub_QuadrOpt -------------------------------------------------
function [stControl, vfNewX, vfNewY, viNewIdx, fCurve] = BSIM_Sub_QuadrOpt(stControl, x, y)
%---------------------------------------------------------------
% function: BSIM_Sub_QuadrOpt
% fits a quadratic equation to the input values and adjusts the points
% according to the estimated maximum
%
% Inputs :  stControl (Structure containing parameters)
%           x (delays of EC-process used for quadratic interpolation)
%           y (SNR values corresponding to the delays)
%           x and y are both of length 3!

% Outputs:  stControl (Structure containing parameters (not changed in this function))
%           vfNewX (delay corresponding to a maximum SNR)
%           vfNewY (estimated maximum SNR)
%           viNewIdx (Idx of the est. Maximum)
%           fcurve (<0: Maximum; >0: Minimum)

% Comments added by Christopher Hauth and Thomas Brand
% Date: 03.Dez.2013
%---------------------------------------------------------------

x=x(:);
y=y(:);

if length(y) ~= 3 || length(x) ~= 3
    error('only three points in quadratic estimation implemented. sorry!');
end

% fit quadratic equation to the three input points
% vfQuadrParm = [(x/stControl.model.EC.fScaleFactor).^2 (x/stControl.model.EC.fScaleFactor) ones(3,1)]\y;
vfQuadrParm = pinv([(x/stControl.model.EC.fScaleFactor).^2 (x/stControl.model.EC.fScaleFactor) ones(3,1)])*y;

% get x/y form from parametric equation
fCurve       = vfQuadrParm(1)/stControl.model.EC.fScaleFactor.^2;
% Estimated delay corresponding to...
fEstExtremeX = (-vfQuadrParm(2)/(2*vfQuadrParm(1)+eps))*stControl.model.EC.fScaleFactor;
% ...the estimated maximum SNR
fEstExtremeY = vfQuadrParm(3)-vfQuadrParm(2)^2/(4*vfQuadrParm(1)+eps);


% we want to find a maximum - so move the point with the smallest y (smallest SNR) close to
% the estimated maximum:
% get index of smallest SNR value
iMinYPos     = find(y == min(y),1,'first');
% These are the other two indices
viNotMinYPos = setdiff(1:3,iMinYPos);
if fCurve < 0     % parabola opens to bottom -> maximum
    % remove the idx of Minimum SNR from delayvector and add idx of est.
    % Maximum SNR
    vfNewX = [fEstExtremeX;x(viNotMinYPos)];
else              % move towards maximum
    if fEstExtremeX < x(1)
        vfNewX  = [fEstExtremeX+stControl.band.EstPer*stControl.model.EC.fSearchRatio;x(viNotMinYPos)];
    elseif fEstExtremeX < x(2)
        viNotMinYPos = [2 3];
        vfNewX  = [fEstExtremeX+stControl.band.EstPer*stControl.model.EC.fSearchRatio;x(viNotMinYPos)];
    elseif fEstExtremeX < x(3)
        viNotMinYPos = [1 2];
        vfNewX  = [fEstExtremeX-stControl.band.EstPer*stControl.model.EC.fSearchRatio;x(viNotMinYPos)];
    else
        vfNewX  = [fEstExtremeX-stControl.band.EstPer*stControl.model.EC.fSearchRatio;x(viNotMinYPos)];
    end
end
% sort x and y values entries by increasing x values
[vfNewX,order] = sort(vfNewX);
vfNewY         = [y(iMinYPos);y(viNotMinYPos)];
vfNewY         = vfNewY(order);
viNewIdx       = find(order == 1);  % only recalculate the new value

% if the test points are too close together, numercial problems could
% arise, so the new point is shifted if necessary
fDivRatio = (vfNewX(2)-vfNewX(1))/(vfNewX(3)-vfNewX(1));
if fDivRatio < stControl.model.EC.fMinPointRatio || fDivRatio > stControl.model.EC.fMaxPointRatio
    if fDivRatio < stControl.model.EC.fMinPointRatio
        r = stControl.model.EC.fMinPointRatio/(1-stControl.model.EC.fMinPointRatio) * (1 + 0.01);%*randn(1));
    else
        r = stControl.model.EC.fMaxPointRatio/(1-stControl.model.EC.fMaxPointRatio) * (1 - 0.01);%*randn(1));
    end
    switch viNewIdx
        case 1
            vfNewX(1) = vfNewX(2) - r * (vfNewX(3)-vfNewX(2));
        case 2
            vfNewX(2) = vfNewX(1) + r/(1+r) * (vfNewX(3)-vfNewX(1));
        case 3
            vfNewX(3) = vfNewX(2) + 1/r * (vfNewX(2)-vfNewX(1));
    end
end


%% BSIM_Sub_CoarseEstMax ---------------------------------------------
function [stControl, fAlpha, fDelay, fMaxSNR] = BSIM_Sub_CoarseEstMax(stControl,spmfXSpec,vfIntensities,spvfOmega)
%---------------------------------------------------------------
% function: BSIM_Sub_CoarseEstMax
% get a coarse estimate of position of the maximum snr by calculating the
% crosscorrelation
% Inputs :  stControl (structure with parameters)
%           spmfXspec (Crosscorrelation Spectra (sparsed matrix))
%           vfIntensities (Intensities of speech and noise in each ear(channel))
%           spvfOmega (sparsed vector containing frequencies (rad) at which the tf of the filter is evaluated)
%
% Outputs : stControl (structure with parameters)
%           fAlpha    (Attenuation, which shall be applied to the left signal due to the EC-Process)
%           fDelay    (Delay, which shall be applied to the left ear signal due to the EC-Process)
%           fMaxSNR   (Maximum SNR, which can be achieved by the EC-process)
% NOTE:
% even though it is said in the article "Revision, extension, and evaluation of a binaural speech
% intelligibility model" by Beutelmann, Brand and Kollmeier (2010) that the tau values are
% restriced to ± 10 ms, this limitation is not included in the code!
% As you will find in the following, the time lag interval depends on the FFT size
% and, more important, on the centre frequency (period) of the gammatone filter 

% Comments added by Christopher Hauth and Thomas Brand
% Date: 03.Dez.2013
%---------------------------------------------------------------
% calculate period of gammatone filter
fPeriod = 1/stControl.model.FB.vfCenterFreqs(stControl.band.iCfIdx);
% limit the timelag according to number of fft points
vfTimeLag  = (-ceil(stControl.signal.nFFT/2-1)-1 : ceil(stControl.signal.nFFT/2-1) ).'/stControl.signal.iFs;
% calculate crosscorrelation
vfXCorr = fftshift(ifft(full(ifftshift(spmfXSpec*stControl.signal.iFs,1))));

% define start and end index to find the maximum of in the cross correlation (within the period of the gammatone filter)
iStartP = find(vfTimeLag >= -fPeriod/2,1,'first');
iEndP   = find(vfTimeLag <=  fPeriod/2,1,'last');
% Find maximum in Crosscorrelation of noise:
[dummy,iNoisePos]  = max(abs(vfXCorr(iStartP:iEndP,2)).*exp(-vfTimeLag(iStartP:iEndP).^2/3e-3));
iNoisePos = iNoisePos + iStartP - 1;
% Again, find the maximum of the Crosscorrelation within the period of the
% gammatone filter, centered around the first estimated maximum
% this time
iStart = find(vfTimeLag >= vfTimeLag(iNoisePos)-fPeriod/2,1,'first');
iEnd   = find(vfTimeLag < vfTimeLag(iNoisePos)+fPeriod/2,1,'last');
[dummy,iNoisePos] = max(real(vfXCorr(iStart:iEnd,2)));
% resulting index for noise:
iNoisePos = iNoisePos + iStart - 1;

% Find maximum in crosscorrelation of speech (same procedure as for noise):
[dummy,iSpeechPos]  = max(abs(vfXCorr(iStartP:iEndP,1)).*exp(-vfTimeLag(iStartP:iEndP).^2/3e-3));
iSpeechPos = iSpeechPos + iStartP - 1;
iStart = find(vfTimeLag >= vfTimeLag(iSpeechPos)-fPeriod/2,1,'first');
iEnd   = find(vfTimeLag < vfTimeLag(iSpeechPos)+fPeriod/2,1,'last');
[dummy,iSpeechPos] = max(real(vfXCorr(iStart:iEnd,1)));
% resulting index for speech
iSpeechPos = iSpeechPos + iStart - 1;

% find the alpha, which maximizes the SNR (tau/delay is given by iSpeechPos/iNoisePos)
% use the time lag found for noise (iNoisePos)
[stControl, vfSpeechLevel, vfNoiseLevel, vfTau, vfAlpha] = BSIM_Sub_ECLevel(stControl,struct('bUseFFT',0,'bOptAlpha',1,'bMonoLimit',0),vfIntensities,[],vfTimeLag(iNoisePos),spvfOmega,spmfXSpec);
% resulting SNR noise
fNoiseSNR = vfSpeechLevel/vfNoiseLevel;
fNoiseAlpha = vfAlpha;

% use the time lag found for speech (iSpeechPos)
[stControl, vfSpeechLevel, vfNoiseLevel, vfTau, vfAlpha] = BSIM_Sub_ECLevel(stControl,struct('bUseFFT',0,'bOptAlpha',1,'bMonoLimit',0),vfIntensities,[],vfTimeLag(iSpeechPos),spvfOmega,spmfXSpec);
% resulting SNR speech
fSpeechSNR = vfSpeechLevel/vfNoiseLevel;
fSpeechAlpha = vfAlpha;
% determine, which procedure (which time lag) leads to the max. SNR  and
% use the corresponding alpha and delay for further calculations
if  fSpeechSNR > fNoiseSNR
    fAlpha  = fSpeechAlpha;
    iDelayPos = iSpeechPos;
    fMaxSNR = fSpeechSNR;
else
    fAlpha  = fNoiseAlpha;
    iDelayPos = iNoisePos;
    fMaxSNR = fNoiseSNR;
end

% add two points around delay position of maximum SNR (three points are needed for quadratic interpolation-> find exact
% location of maximum)
switch iDelayPos
    case 1
        iDelayPos = [1 2 3];
    case length(vfTau)
        iDelayPos = length(vfTau) - [2 1 0];
    otherwise
        iDelayPos = iDelayPos + [-1 0 1];
end
fDelay  = reshape(vfTimeLag(iDelayPos),1,[]);


%% BSIM_Sub_ECLevel --------------------------------------------
function [stControl, vfSpeechLevel, vfNoiseLevel, vfTau, vfAlpha] = BSIM_Sub_ECLevel(stControl,stFlags,vfIntensities,vfAlpha,vfTau,spvfOmega,spmfXSpec)
%---------------------------------------------------------------
% function: BSIM_Sub_ECLevel 
% bFlags: bUseFFT bOptAlpha bMonoLimit
% Realisation of Equilization Cancellation Process
% Inputs:    stControl (structure containing parameters)
%            stFlags (structur containing flags)
%            vfIntensities (Intensities)
%            vfAlpha (Attenuation, which shall be applied in the EC process)
%            (needs to be calculated in this function)
%            vfTau (Delay, which shall be appied in the EC process)
%            spvfOmega (sparsed vector containing frequencies (rad) at which the tf of the filter is evaluated)
%            spmfXSpec (sparsed matrix containing Crosscorrelation Spectra)

% Outputs:  stControl
%           vfSpeechLevel (Level of Speech after EC-processing)
%           vfNoiseLevel  (Level of Noise after EC-processing)
%           vfTau         (Delay, which has been applied in the EC-process; Output for further processing)
%           vfAlpha       (Attenuation, which has been applied in "")

% Comments added by Christopher Hauth and Thomas Brand
% Date: 03.Dez.2013
%--------------------------------------------------------------------------

if ~isfield(stFlags,'bUseFFT')
    stFlags.bUseFFT = false;
end
if ~isfield(stFlags,'bOptAlpha')
    stFlags.bOptAlpha = false;
end
if ~isfield(stFlags,'bMonoLimit')
    stFlags.bMonoLimit = false;
end

if stFlags.bUseFFT
    vfTimeLag  = ( -ceil(stControl.signal.nFFT/2-1)-1 : ceil(stControl.signal.nFFT/2-1) ).'/stControl.signal.iFs;
    [dummy,iTimeStart] = min(abs(vfTimeLag-min(vfTau)));
    [dummy,iTimeEnd  ] = min(abs(vfTimeLag-max(vfTau)));
    
    vfTau   = vfTimeLag(iTimeStart:iTimeEnd);
    
    vfXCorr = fftshift(- 2 .* real(ifft(full(ifftshift(spmfXSpec*stControl.signal.iFs,1)))),1);
    vfXCorr = vfXCorr(iTimeStart:iTimeEnd,:);
    iNumPoints = length(vfTau);
else
    iNumPoints = length(vfTau);
    viNZIdx = find(spvfOmega); % get indices out of sparsed matrix
    % calculate standard deviation of delay error according to equation 1 in
    % Beutelmann R, Brand T (2006), section 4: Artificial processing errors
    vfSigmaDelta = stControl.model.EC.fSigmaDelta0 * ( 1 + abs(vfTau)/stControl.model.EC.fDelta0 );
    if ~stControl.model.EC.bErrOn % if no processing errors are used, standard deviation is set to zero
        vfSigmaDelta(:) = 0;
    end
    % sets a lower limit for correlation decrement at high frequencies
    % additional factor of [log(10)/20] results from dB(basis 10) to neper scale (basis e)
    % conversion (errors are used in exponential function)
    fLowLimit = stControl.model.EC.fLowLimit0*exp((log(10)/20*stControl.model.EC.fSigmaEpsilon0).^2);
    if ~stControl.model.EC.bErrOn
        fLowLimit(:) = 0;
    end
    % for the following equations see Appendix A: Detailed Derivation (Beutelmann, R., Brand, T. and Kollmeier, B. (2010))
    % Computation of overall intensity of the residual speech and noise
    % 1) build matrix z containing delay errors (see cross correlation term of equation A5 )
    %   exp(-spvfOmega(viNZIdx).^2*vfSigmaDelta.^2) = gaussian lowpass filter applied on cross correlation term
z = sparse(repmat(viNZIdx,1,iNumPoints),repmat(1:iNumPoints,length(viNZIdx),1),...
               ((1-fLowLimit)*exp(-spvfOmega(viNZIdx).^2*vfSigmaDelta.^2)+fLowLimit)...% Gaussian lowpass
               .* exp(1i*spvfOmega(viNZIdx)*vfTau),...% phase (xy*)
               stControl.signal.nFFT,iNumPoints);

% 2) Compute cross correlation term  (see cross correlation term of equation A5, Beutelmann, R., Brand, T. and Kollmeier, B. (2010)) 
for ch = 1:2                                                                 
        vfXCorr(:,ch) = -2 * real(sum( repmat(spmfXSpec(:,ch),1,length(vfTau)).* z)).' * stControl.signal.iFs/stControl.signal.nFFT;
end

end

if stFlags.bOptAlpha
    bOptAlpha = true;
    
% 3) calculate optimal alphas: analytical solution of equation A5:
   % for detailed information about the derivation of the following
   % equations, see documentation of BSIM
    fSqEqSolvePart1 = vfIntensities(2) * vfIntensities(3) - vfIntensities(1)  * vfIntensities(4);

    fSqEqSolvePart2 = vfIntensities(1) * vfXCorr(:,2) - vfIntensities(3) * vfXCorr(:,1) + eps;
    fSqEqSolvePart3 = vfIntensities(2) * vfXCorr(:,2) - vfIntensities(4) * vfXCorr(:,1) + eps;

    % calculate alpha analytically  by using the equations above in the pq-formular:
    vfAlpha = (fSqEqSolvePart1 + sqrt(fSqEqSolvePart1^2+fSqEqSolvePart2.*fSqEqSolvePart3)*[-1 1])./repmat(fSqEqSolvePart2,1,2);
    vfAlpha(vfAlpha <= 0) = 1;
end
% standard deviation of gain error according to equation (1), Beutelmann R, Brand T (2006): 
% Deviating from the formular, a factor of (log10/20) is used to transform the dB scaling into 
% neper scaling -> ( basis e instead of 10)
vfSigmaEpsilon = (log(10)/20) * stControl.model.EC.fSigmaEpsilon0 * ( 1 + ( abs(20*log10(vfAlpha)) / stControl.model.EC.fAlpha0 ).^stControl.model.EC.fp );
if ~stControl.model.EC.bErrOn
    vfSigmaEpsilon(:) = 0; % standard deviation of gain error is set to zero, if no processing errors are used
end
% apply alpha to left ear channel:
% exp(-vfSigmaEpsilon.^2) is the gain error
vfSpeechLevel = ( vfAlpha.^2 * vfIntensities(1) + vfIntensities(2) + vfAlpha .* exp(-vfSigmaEpsilon.^2) .* repmat(vfXCorr(:,1),1,size(vfAlpha,2))) ./ ((vfAlpha > 1) .* vfAlpha.^2 + (vfAlpha <= 1));
vfNoiseLevel  = ( vfAlpha.^2 * vfIntensities(3) + vfIntensities(4) + vfAlpha .* exp(-vfSigmaEpsilon.^2) .* repmat(vfXCorr(:,2),1,size(vfAlpha,2))) ./ ((vfAlpha > 1) .* vfAlpha.^2 + (vfAlpha <= 1));

if stFlags.bOptAlpha
    % select best snr
    mfSNR = vfSpeechLevel./vfNoiseLevel;
    [mfSNR,miMaxIdx] = max(mfSNR,[],2);
    vfAlpha = vfAlpha([miMaxIdx==1 miMaxIdx==2]);
    vfSpeechLevel  =  vfSpeechLevel([miMaxIdx==1 miMaxIdx==2]);
    vfNoiseLevel   =   vfNoiseLevel([miMaxIdx==1 miMaxIdx==2]);
    vfSigmaEpsilon = vfSigmaEpsilon([miMaxIdx==1 miMaxIdx==2]);
end
% monaural SNR (better ear listening)
[fMaxMonoSNR,iMaxMonoIdx] = max(vfIntensities([1 2])./vfIntensities([3 4]));
vfSpeechLevel(isnan(vfSpeechLevel)) = vfIntensities(iMaxMonoIdx);
vfNoiseLevel(isnan(vfSpeechLevel))  = vfIntensities(iMaxMonoIdx+2);

if stFlags.bMonoLimit
    % limit all values to be greater than the best monaural snr
    iLowIdx = vfSpeechLevel./vfNoiseLevel < fMaxMonoSNR;
    vfSpeechLevel(iLowIdx) = vfIntensities(iMaxMonoIdx)   * (1-0.01*abs(vfTau(iLowIdx))/max(abs(vfTau)+eps));
    vfNoiseLevel(iLowIdx)  = vfIntensities(iMaxMonoIdx+2);
end

if isfield(stControl.band,'bSaveNow') && stControl.band.bSaveNow
    stControl.band.vfInt   = vfIntensities(:);
    stControl.band.vfCross = full(exp(-vfSigmaEpsilon.^2) .* vfXCorr);
    stControl.band.bSaveNow = false;
end

%% BSIM_Sub_TransferFunction -----------------------------------------
function [stControl, spvfTransferFct, spvfOmega] = BSIM_Sub_TransferFunction(stControl,fBwRatio)
%--------------------------------------------------------------------------
% function: BSIM_Sub_TransferFunction
% depending on the realisation of the filterbank mirroring the frequency
% selectivity on the basilar membrane,
% filter function subroutines are used
% Input:    stControl (Structure containing parameters)
%           fBwRatio  (Bandwidth ratio = 1)

% Output:   stControl (Structure containing parameters)
%           spvfTransferFct (sparsed vector cont. TF of Filter)
%           spvfOmega (sparsed vector cont. frequencies (rad) at which TF is evaluated)

% Comments added by Christopher Hauth and Thomas Brand
% Date: 03.Dez.2013
%--------------------------------------------------------------------------

switch lower(stControl.model.sFilterType)
    case 'gtf' % gammatone filter
        [stControl, spvfTransferFct, spvfOmega] = BSIM_Sub_GammatoneFilter(stControl,fBwRatio);
    case 'third' % third-octave filter
        [stControl, spvfTransferFct, spvfOmega] = BSIM_Sub_ThirdOctaveFilter(stControl);
    case 'cb' % critical band
        [stControl, spvfTransferFct, spvfOmega] = BSIM_Sub_CriticalBandFilter(stControl);
    otherwise
        error('unknown filter type. filter types are: gtf (gammatone filter), third (third octave bandpass), cb (critical band bandpass)');
end



%% BSIM_Sub_CriticalBandFilter ----------------------------------------
function [stControl, spvfTransferFct, spvfOmega] = BSIM_Sub_CriticalBandFilter(stControl)
%---------------------------------------------------------------
% function: BSIM_Sub_CriticalBandFilter
% function to generate a Filter with critical bandwidth

% Input:    stControl (Structure containing parameters)

% Output:   stControl (Structure containing parameters)
%           spvfTransferFct (sparsed vector cont. TF of Critical Band Filter)
%           spvfOmega (sparsed vector cont. frequencies (rad) at which TF is evaluated)

% Comments added by Christopher Hauth and Thomas Brand
% Date: 03.Dez.2013 
%---------------------------------------------------------------

% bark scale (from ANSI S3.5-1997 (SII))
vfCf   = [150 250 350 450 570 700 840 1000 1170 1370 1600 1850 2150 2500 2900 3400 4000 4800 5800 7000 8500];
vfLow  = [100 200 300 400 510 630 770  920 1080 1270 1480 1720 2000 2320 2700 3150 3700 4400 5300 6400 7700];
vfHigh = [200 300 400 510 630 770 920 1080 1270 1480 1720 2000 2320 2700 3150 3700 4400 5300 6400 7700 9500];

% frequency scale
vfFreq = ((-floor(stControl.signal.nFFT/2):ceil(stControl.signal.nFFT/2)-1)/stControl.signal.nFFT).' * stControl.signal.iFs;
viLimits = find(vfFreq >= vfLow(stControl.band.iCfIdx) & vfFreq < vfHigh(stControl.band.iCfIdx));
spvfOmega = sparse(viLimits,1,2*pi*vfFreq(viLimits),stControl.signal.nFFT,1);

% transfer function
spvfTransferFct = spones(spvfOmega);

% filter bandwidth in hz
stControl.model.FB.vfBandWidth(stControl.band.iCfIdx,1) = vfHigh(stControl.band.iCfIdx)-vfLow(stControl.band.iCfIdx);


%% BSIM_Sub_ThirdOctaveFilter ----------------------------------------
function [stControl, spvfTransferFct, spvfOmega] = BSIM_Sub_ThirdOctaveFilter(stControl)
%---------------------------------------------------------------
% function: BSIM_Sub_ThirdOctaveFilter
% function to create a 1/3-Octave filter
% 
% Input:    stControl (Structure containing parameters)


% Output:   stControl (Structure containing parameters)
%           spvfTransferFct (sparsed vector cont. TF of 1/3-octave Filter)
%           spvfOmega (sparsed vector cont. frequencies (rad) at which TF is evaluated)

% Comments added by Christopher Hauth and Thomas Brand
% Date: 03.Dez.2013
%---------------------------------------------------------------

% frequency scale
vfFreq = ((-floor(stControl.signal.nFFT/2):ceil(stControl.signal.nFFT/2)-1)/stControl.signal.nFFT).' * stControl.signal.iFs;
viLimits = find(vfFreq >= stControl.model.FB.vfCenterFreqs(stControl.band.iCfIdx)*2^(-1/6) & vfFreq < stControl.model.FB.vfCenterFreqs(stControl.band.iCfIdx)*2^(1/6));
spvfOmega = sparse(viLimits,1,2*pi*vfFreq(viLimits),stControl.signal.nFFT,1);

% transfer function
spvfTransferFct = spones(spvfOmega);

% filter bandwidth in hz
stControl.model.FB.vfBandWidth(stControl.band.iCfIdx,1) = stControl.model.FB.vfCenterFreqs(stControl.band.iCfIdx)*2^(1/3);



%% BSIM_Sub_GammatoneFilter ------------------------------------------
function [stControl, spvfTransferFct, spvfOmega] = BSIM_Sub_GammatoneFilter(stControl,fBwRatio)
%---------------------------------------------------------------
% function: BSIM_Sub_GammatoneFilter
% function to generate a Gammatonefilter
% 
% Input:    stControl (Structure containing parameters)
%           fBwRatio  (Bandwidth ratio = 1)

% Output:   stControl (Structure containing parameters)
%           spvfTransferFct (sparsed vector cont. TF of Gammatonefilter)
%           spvfOmega (sparsed vector cont. frequencies (rad) at which TF is evaluated)

% Comments added by Christopher Hauth and Thomas Brand
% Date: 03.Dez.2013
%---------------------------------------------------------------
% construct a gammatone filter of width fBwRatio/ERB at stControl.model.FB.vfCenterFreqs(stControl.band.iCfIdx)/Hz:
% Realisation of gammatone filter according to  Hohmann, 2002, Acta Acoustica 88


% calculate bandwidth in ERB at stControl.model.FB.vfCenterFreqs(stControl.band.iCfIdx)
fBw_ERB = fBwRatio * (stControl.model.FB.fERBOffset + ...
    stControl.model.FB.vfCenterFreqs(stControl.band.iCfIdx)/stControl.model.FB.fERBScale);
% bandwidth factor
a = (pi*factorial(2*stControl.model.FB.iFilterOrder-2)*2^(-(2*stControl.model.FB.iFilterOrder-2)))/(factorial(stControl.model.FB.iFilterOrder-1)^2);
b = fBw_ERB/a;
fLambda = exp(-2*pi*b/stControl.signal.iFs);
% centre frequency factor
fBeta = 2*pi*stControl.model.FB.vfCenterFreqs(stControl.band.iCfIdx)/stControl.signal.iFs;
% filter coefficient
fATilde = fLambda * exp(1i*fBeta);
% normalization factor
fNormFactor = 2 * (1 - abs(fATilde)) ^ stControl.model.FB.iFilterOrder;

% find indices between which the transfer function will be greater than
% stControl.GTB.fZeroThresh, i.e. significantly above from zero (for out
% purposes)
viLimits = floor(stControl.signal.nFFT/2) + 1 ...   %
    + round( (fBeta + [-1 1] * acos(fLambda/2 + 1/(2*fLambda) ...
              * ( 1 - (10^(stControl.model.FB.fZeroThresh/20)/fNormFactor)^(-2/stControl.model.FB.iFilterOrder) )) ...
             )*stControl.signal.nFFT/2/pi);
viLimits(1) = max([1,viLimits(1)]);        % in case something goes wrong numerically
viLimits(2) = min([stControl.signal.nFFT,viLimits(2)]);  %                - " -

% frequency scale
% full scale would be: 2*pi((-floor(stControl.signal.nFFT/2):ceil(stControl.signal.nFFT/2)-1)/stControl.signal.nFFT).';
spvfOmega = sparse(viLimits(1):viLimits(2),1,2*pi*(-floor(stControl.signal.nFFT/2)-1+(viLimits(1):viLimits(2)))/stControl.signal.nFFT,stControl.signal.nFFT,1);

% transfer function
% factor of 0.5 in order to maintain the correct level
spvfTransferFct = spfun(@(om)0.5*fNormFactor./((1-fATilde*exp(-1i*om)).^stControl.model.FB.iFilterOrder),spvfOmega);
% convert omega to actual units for later use
spvfOmega = spvfOmega * stControl.signal.iFs;

% filter bandwidth in Hz
stControl.model.FB.vfBandWidth(stControl.band.iCfIdx,1) = 2*sqrt(2^(1/stControl.model.FB.iFilterOrder)-1)*b;



%% BSIM_Sub_Constants ------------------------------------------------
function BSIM_Sub_Constants
%---------------------------------------------------------------
% function: BSIM_Sub_Constants
% creates constants in the main function's workspace

% Comments added by Christopher Hauth and Thomas Brand
% Date: 03.Dez.2013
%---------------------------------------------------------------

% Add simple constants to this cell array (in the format shown in the remark).
% They will automatically be created in the main function's workspace.
caConstants = {
    %     'constant name','value','extoverride';...
    'stControl.run.sDisplay',      '''text''',          1; ... % display information type
    'stControl.signal.iFs',         '44100',            1; ...
    'stControl.HL.mfAudiogram', '[0 0 0;50000 0 0]',    1; ... % standard audiogram , 0 dB HL for all frequencies
    'stControl.HL.fdBDetThresh',        '-1',          1; ...  % (~1?) dB, threshold criterion for the internal noise
    'stControl.model.sFilterType',  '''gtf''',          1; ... % type of filter bank: gammatone (gtf), third octave (third)
    'stControl.model.FB.iFilterOrder',  '4',            1; ... % gammatone fb order (4th order)
    'stControl.model.FB.fERBOffset',   '24.7',          1; ... % erb scale constants (Hohmann, 2002, acta acoustica 88)
    'stControl.model.FB.fERBScale',     '9.265',        1; ... %      - " -                     - " -
    'stControl.model.FB.fZeroThresh',   '-20',          1; ... % (1.139,2.288562) dB, disregard parts of the transfer fct which are lower
    'stControl.model.FB.fERBBwRatio',   '1',            1; ... % standard bandwidth ratio
    'stControl.model.EC.fMaxDelay',     '1e-2',         1; ... % (0.003492) sec., maximum delay (for zero padding)
    'stControl.model.EC.bErrOn',        'false',        1; ... % controls processing errors of EC mechanism (Beutelmann R, Brand T (2006)):
                                                               % parameters controlling standard deviation of delay errors: see formula(1) 
    'stControl.model.EC.Instant_flank', '1',            1; ... % Instant BSIM flank, default: 1
    'stControl.model.EC.fSigmaDelta0',  '65e-6',        1; ... % default: 65e-6 [s]
    'stControl.model.EC.fDelta0',       '1.6e-3',       1; ... % default: 1.6e-3 [s] 
                                                               % parameters controlling standard deviation of gain errors: see formula(1)
    'stControl.model.EC.fSigmaEpsilon0','1.5',          1; ... % default: 1.5 [dB]
    'stControl.model.EC.fAlpha0',      '13.0',          1; ... % default: 13 [dB]
    'stControl.model.EC.fp',            '1.6',          1; ... % default: 1.6 (no unit->used as exponential factor)
                                                               % all parameters were calculated by vom H�vel(1984) 
    
    'stControl.model.EC.fLowLimit0',    '1/3',          1; ... % 1/3 linear lower limit of correlation low pass for high frequencies
    'stControl.model.EC.fScaleFactor',  '1.0e-3',       1; ... % scale factor can be used to avoid numerical problems in the optimum search algorithm
    'stControl.model.EC.fSearchRatio', '(sqrt(5)-1)/2', 1; ... % ratio for distribution of points
	'stControl.model.EC.fMinPointRatio','0.05',         1; ...
    'stControl.model.EC.fMaxPointRatio','0.95',         1; ...
    'stControl.model.LastPar.bIsSet'   ,'false',        1; ... % binaural sluggishness: last value is valid
    };

csFilterTypes = {'gtf','third','cb'};
% gammatone filterbank centre frequencies
csCenterFreqs{1} = '[146.02 188.74 236.34 289.36 348.42 414.21 487.50 569.15 660.10 761.41 874.27 1000.00 1140.06 1296.07 1469.87 1663.48 1879.16 2119.41 2387.05 2685.19 3017.31 3387.29 3799.43 4258.55 4769.99 5339.72 5974.39 6681.39 7468.98 8346.32].'''; 
% one-third octave band SII procedure (ANSI S3.5-1997 Table 3)
csCenterFreqs{2} = '[    160    200     250   315    400            500    630           800      1000         1250            1600            2000     2500                  3150             4000            5000            6300            8000].''';       
% critical band SII procedure (ANSI S3.5-1997 Table 1)
csCenterFreqs{3} = '[  150              250         350      450       570        700     840     1000     1170       1370     1600      1850      2150 2500            2900       3400        4000         4800         5800          7000          8500].'''; 

iFilterTypeIdx = find(strcmp(csFilterTypes,caConstants{find(strcmp(caConstants(:,1),'stControl.model.sFilterType')),2}(2:end-1)));
caConstants(end+1,:) = {'stControl.model.FB.vfCenterFreqs',csCenterFreqs{iFilterTypeIdx},1}; % standard (gammatone,third,critical band) filter bank center frequencies (140-9000 Hz, 1000Hz+/-xERB, 160-8000 Hz, 150-8500 Hz)

% create constants in main function workspace
for k = 1:size(caConstants,1)
%     disp(caConstants{k,1});
    if caConstants{k,3}
        try
            val = evalin('caller',caConstants{k,1}); % create predefined variable here
        catch
            evalin('caller',sprintf('%s = %s;',caConstants{k,1},caConstants{k,2})); % set to default otherwise
        end
    else
%         assignin('caller',caConstants{k,1},caConstants{k,2});
        evalin('caller',sprintf('%s = %s;',caConstants{k,1},caConstants{k,2}));
    end
end

%% Add more complex constant definitions here:
%% Include hearing threshold as defined in Table 1 (Tf) ISO 226, 2003
%% (%compute dB SPL from dB HL (ISO 226, 2003)-RainerBeutelmann)

% changed reference to match approximately SII reference internal noise
% spectrum, see also internal noise calc procedure (RainerBeutelmann):
% Values were changed according to ANSI Methods for Calculation of the
% Speech Intelligibility Index (June 1997), p.7, section 3.24
% Reference internal noise spectrum level, Note 2:
% the internal noise spectrum level is equal to: 

% Threshold = pure tone threshold (as defined in ISO 226, 2003) - 10log10(ERBn/1Hz),

% where ERBn = 24.7*(4.37*F(kHz)+1);  see: Glasberg and Moore, 1990
%       Equivalent Rectangular Bandwidth
% resulting values of the changed reference:
evalin('caller',sprintf('%s = %s;','stControl.HL.mfFreqHL2SPL',...
 '[20,62.41;25,52.92;32,43.91;40,35.47;50,28.41;63,21.71;80,15.47;100,10.20;125,5.48;160,0.87;200,-2.96;250,-6.44;315,-9.69;400,-12.62;500,-15.16;630,-17.27;800,-18.66;1000,-19.03;1250,-19.13;1600,-21.46;2000,-25.01;2500,-28.89;3150,-31.92;4000,-32.70;5000,-30.02;6300,-22.98;8000,-15.99;10000,-15.63;12500,-15.48]'));


% ******************************************************************************
% * helper functions ***********************************************************
% ******************************************************************************

%% dB ----------------------------------------------------------------
function fOut_dB = dB(fIn_linear)
%---------------------------------------------------------------
% function: dB
% calculates dB from linear values 
%---------------------------------------------------------------
fOut_dB = 20*log10(abs(fIn_linear));

% end of file