function output_signal = simulate_hearing_loss(signal,fs,hearing_thresholds,doPlot,pegel_cor)
%simulate_hearing_loss - simulate HI loudness perception 
%purpose: 
%
% Syntax:  output_signal = simulate_hearing_loss(signal,fs,hearing_thresholds)
% Inputs:
%   signal             - input mono signal
%   fs                 - sampling frequency, should be 44,100 Hz, others 
%                        might work but it was not carefully checked
%   hearing_thresholds - hearing thresholds in dB HL at 250, 500, 1000,
%                        2000, 4000, 8000 Hz
%   doPlot             - set to 1 to show debeug plots
%   pegel_cor          - Variable to adjustment the Pegel 1-3 ISO 226, 4-6 
%                        iso226:2003(E), 7 iso389-8:2003 Table C.1 (default 4)
%
% Outputs:
%   output_signal - signal modified based on the monaural loudness functions
%
%
% Example:
%    Line 1 of example
%    Line 2 of example
%
% -------------------------------------------------------------------------
%                    2016 Fraunhofer IDMT, Oldenburg   (Copyright)       
%
% This software and/or program is protected by copyright law and international
% treaties. Any reproduction or distribution of this software and/or program,
% or any portion of it, may result in severe civil and criminal penalties,
% and will be prosecuted to the maximum extent possible under law.
% -------------------------------------------------------------------------
% author          : Dirk Oetting          
% email           : Dirk.Oetting@idmt.fraunhofer.de           
% version         : 0.1         
% date            : 11-May-2016           
% Matlab          : 7.9.0.529 (R2009b)   
% HSA style guide : 0.3
% -------------------------------------------------------------------------
% changes  : -    11.5.2016 Implementation of basic functionality
%            -    13.5.2016 added windowing of the output signal (DO)
%            -    26.5.2016 fixed fft/ifft synthesis. Now the error is in
%            the domain of 1e-11 dB, if zero gain is applied.
%            -    21.06.2016 Variable to adjustment the Pegel in different
%                 ways
%
if nargin < 5
    pegel_cor = 4;
elseif nargin < 4 
pegel_cor = 4;
doPlot = 0;
end
    % hearing_thresholds in dB HL at 250, 500, 1000, 2000, 4000, 8000 Hz
    if size(signal,2)~=1
        warning('need mono signal, only taking first channel');
        signal = signal(:,1);
    end
    
    % extend signal length to an even value
    if rem(length(signal),2)
        signal(end+1) = 0;
    end
    
    if fs~=44100
        warning('sample rate should be 44,100 Hz');
    end

    samplerate = 44100;				% sampling rate in Hz
    center_f = [250 500 1000 2000 4000 8000];    
    
    % Reference loudness function for narrow-band signals
    % in dB HL
    % frequency, Lcut (dB HL), m_low, m_high
    % fitparams chosen that HTL = 0dB and UCL = 100 dB
    % with Lcut at 85 dB and mlow around 0.25. Compare
    % data of Oetting et al. (2015) Spectral and binaural 
    % loudness summation for hearing-impaired listeners, Hear Res 335
    % 
    % mean data of Oetting et al. (2016)
    mean_UCL = 107.9; % dB HL    
    p_fit_HTL_Lcut = [0.2843   77.0525];
    fitparamsNH = [];
    for idxFreq = 1:length(center_f)
        HTL = 0;
        mlow = 0.29;
        CP = 25;
        b = 2.5 - mlow*HTL;
        Lcut = (CP - b)/mlow;
        UCL = mean_UCL;
        mhigh = (50-CP)/(UCL - Lcut);
        fitparamsNH(idxFreq,:) = [Lcut mlow mhigh];
    end

    fitparamsHI = [];
    for idxFreq = 1:length(center_f)
        HTL = hearing_thresholds(idxFreq);
        % relationship between HTL and mlow
        Lcut =  polyval(p_fit_HTL_Lcut,HTL);
        mlow = 22.5/(Lcut-HTL);
        UCL = mean_UCL;
        if Lcut - 5 >= UCL 
            disp('warning: wrong setting for UCL, function converted to linear');
            mhigh = 5;
        else
            mhigh = (50-CP)/(UCL - Lcut);
        end

        fitparamsHI(idxFreq,:) = [Lcut mlow mhigh];
    end

    signal = signal(:);
    % make sure the signal length is even
    if mod(length(signal),2)==1
        signal = [signal 0];
    end
    
    % amplify signal according to NB loudness scaling data
    input_spectrum  = fft(signal)/length(signal);
    
    %Just keep right sided spectrum, e.g. left half in matlab implementation
    input_spectrum = input_spectrum(1:length(signal)/2+1); 
    % spectrum must be multiplied squared bins by two to keep the total energy
    input_spectrum = input_spectrum *  2^(1/2);

    frequency_vector = (0:length(signal)/2) /length(signal) *  samplerate ;
    % set dc bin  bin to zero
    input_spectrum_reduced = [0;input_spectrum(2:end)];
    

    % frequency channels
    edge_f = 10.^(mean([log10(center_f(1:end-1));log10(center_f(2:end))]));
    edge_f = [0 edge_f max(frequency_vector)];
    channel_level_dB_SPL_FF = zeros(length(center_f),1);
    channel_level_dB_HL  = zeros(length(center_f),1);
    channel_gain_dB  = zeros(length(center_f),1);
    CU  = zeros(length(center_f),1);
    
    for idxChannel = 1:length(center_f)
        channel_level_dB_SPL_FF(idxChannel) = ... 
        10*log10(sum(abs(input_spectrum_reduced(frequency_vector<=edge_f(idxChannel+1) & frequency_vector>edge_f(idxChannel))).^2));
        channel_level_dB_HL(idxChannel) = channel_level_dB_SPL_FF(idxChannel) - isothr_oetting(center_f(idxChannel),pegel_cor);    

        %%%CU(idxChannel) = loudness_function_bh2002(channel_level_dB_HL(idxChannel),fitparamsNH(idxChannel,:));
        CU(idxChannel) = loudness_function_bh2002(channel_level_dB_HL(idxChannel),fitparamsHI(idxChannel,:));
        % CU will be between 0 and 50 CU. look into NH and HI loudness function 
        % which gain is needed
        channel_gain_dB(idxChannel) = loudness_function_bh2002(CU(idxChannel),fitparamsHI(idxChannel,:),true) - loudness_function_bh2002(CU(idxChannel),fitparamsNH(idxChannel,:),true);
    end

    gain_dB = interp1(log([0 center_f max(frequency_vector)]+eps),channel_gain_dB([1 1:end end],1),log(frequency_vector+eps))';

    % inverse gain to simulate HI for NH listeners
    gain_dB = - gain_dB;
    % apply gain
    output_spectrum = 10.^(gain_dB/20) .* input_spectrum;	

    % divide by sqrt(2) to keep energy
    output_spectrum  = output_spectrum / 2^(1/2);				

    % That's all what is left to do for calculating the output spectrum
    output_spectrum = [output_spectrum; conj(output_spectrum(end-1:-1:2))];
    output_signal = ifft(output_spectrum)   * length(signal);
    if any(imag(output_signal))
        output_signal = real(output_signal);
        warning('output signal contained imaginary parts');
    end
    
    
    
    
    %% ###################### debug plots ##########################
    if doPlot
        clf
        p = [];
        set(gcf,'WindowStyle','Docked');
        p(1) = plot(frequency_vector,20*log10(abs(input_spectrum)),'color','k');
        set(gca,'XScale','log')
        f_scale = [50 100 250 500 1000 2000 4000 8000 16000];
        set(gca,'xtick',f_scale)
        set(gca,'xticklabel',f_scale)
        hold on;
        for idxChannel = 1:length(center_f)
            plot(edge_f(idxChannel+1)*[1 1],[-200 200],'k');
            p(2) = plot([edge_f(idxChannel)+eps edge_f(idxChannel+1)],[1 1]*channel_level_dB_SPL_FF(idxChannel),'-b','linewidth',4);
            x1 = [1.2 1.2 1.2 1 1 1]*center_f(idxChannel);
            x2 = 1./[1.2 1.2 1.2 1 1 1]*center_f(idxChannel);

            y_NH = loudness_function_bh2002([2.5 25 50],fitparamsNH(idxChannel,:),true) + isothr_oetting(center_f(idxChannel),4);y_NH = [y_NH fliplr(y_NH)];
            y_HI = loudness_function_bh2002([2.5 25 50],fitparamsHI(idxChannel,:),true)+ isothr_oetting(center_f(idxChannel),4);y_HI = [y_HI fliplr(y_HI)];
            f = fill(x1,y_NH,'b');
            set(findall(f,'FaceColor','b'),'FaceColor','y');
            p(3) = f;
            f = fill(x2,y_HI,'r');
            set(findall(f,'FaceColor','r'),'FaceColor','g');
            p(4) = f;
             plot(x1(3:4),y_NH(2)*[1 1],'-k');
            plot(x2(3:4),y_HI(2)*[1 1],'-k');               

        end
        plot(center_f,isothr_oetting(center_f,4),'or','markersize',12);
        plot(center_f,channel_gain_dB,'xr','markersize',12);

        p(5) = plot(frequency_vector,isothr_oetting(frequency_vector,4),'-g');
        p(6) = plot(frequency_vector,gain_dB,'-r');
        plot(frequency_vector,gain_dB*0,'--r');

        legend(p,'Spektrum','Channel Level','NH CLS','SH CLS','HTL','Gain','Location','NorthWest');
        ylabel('Level / dB SPL FF')
        xlabel('Frequency / Hz')
        grid off;
        xlim([20 18000])
        ylim([-30 120])
    end   











    
    
    

