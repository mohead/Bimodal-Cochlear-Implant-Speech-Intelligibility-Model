function analyse_in_out( sound_one, sound_two, fs, name )
% function analyse_in_out( sound_one, sound_two, fs, name )
% This function plots an analysis of 'sound_one' vs. 'sound_two' with a given
% samplingrate 'fs' in a figure with a title 'name'. More specifically,
% you will see the cross-correlation between the signals, their
% differences in the time-domain, the signals over time + global RMS in dB
% FS, the Power spectrum density (using pwelch) and the spectrogramms for
% both signals.

if length(sound_one) ~= length(sound_two)
[sound_one, sound_two] = append_zeros_to_shorter_signal(sound_one, sound_two); 
warning('Input signals have different lengths, zero-padding last signal')
end

vTime = linspace(0,length(sound_one)/fs,length(sound_one));

vDiff = sound_one - sound_two;

h = figure('Name',name,'Units','normalized','Position',[0.2 0.2 0.6 0.6]); 
   
subplot(2,3,1);
plot(vTime,vDiff)
title('Difference soundone - soundtwo')
xlabel('Time(sek.)')
ylabel('Value')

vCorr = xcorr(sound_one, sound_two);
corr_leg = (2*length(sound_one))/2;

vTimeCorr = [ linspace(-length(sound_one)/fs,0,length(sound_one)) vTime(2:end)];

subplot(2,3,2);
plot(vTimeCorr,vCorr)
title('Crosscorrelation')
xlabel('Time(sek.)')
ylabel('Value')

rms_one = 10*log10(rms2(sound_one,1)+eps);
rms_two = 10*log10(rms2(sound_two,1)+eps);

legen1=sprintf('Sound one %0.2fdBFs',rms_one);
legen2=sprintf('Sound two %0.3fdBFs',rms_two);

subplot(2,3,3)
if rms_one >= rms_two
plot(vTime,sound_one,'-b')
hold on
plot(vTime,sound_two,'-.r')
legend(legen1,legen2)
else
    plot(vTime,sound_two,'-.r');
    hold on;
    plot(vTime,sound_one,'-b');
    legend(['Sound two ' num2str(rms_two) 'dbFs'],['Sound one ' num2str(rms_one) 'dbFS'])
end
title('Signal + RMS')
xlabel('Time(sek.)')


% playback-buttons
uicontrol('Parent',h,'Style', 'pushbutton',...
       'String', 'Play Sound 1',... %replace with the text you want
       'Units','normalized',...
       'Position', [0.91 0.8 0.08 0.08],...
       'Callback', {@callback_play,sound_one,fs});
   
uicontrol('Parent',h,'Style', 'pushbutton',...
       'String', 'Play Sound 2',... %replace with the text you want
       'Units','normalized',...
       'Position', [0.91 0.7 0.08 0.08],...
       'Callback', {@callback_play,sound_two,fs});
   
NFFT = 4096;
[pxx,f] = pwelch(sound_one,NFFT,0,NFFT,fs);
[pyy,f] = pwelch(sound_two,NFFT,0,NFFT,fs);

rms_pwelch(1) = 10*log10(sqrt(sum(pxx)));
rms_pwelch(2) = 10*log10(sqrt(sum(pyy)));

legen1=sprintf('Sound one %0.2fdBFs',rms_pwelch(1));
legen2=sprintf('Sound two %0.3fdBFs',rms_pwelch(2));
subplot(2,3,4);
cor = 5.6992; % max(20*log10(pxx)) - rms_one;
plot(f,10*log10(pxx))
hold on
plot(f,10*log10(pyy))
legend(legen1,legen2)
title('pWelch for both signals')
xlabel('Frequency')
ylabel('Magnitude (dB)')
set(gca,'xscale','log'); %Logarithmic frequency axis

diff = pxx - rms_one;

subplot(2,3,5)
spectrogram(sound_one,NFFT,0,NFFT,fs,'yaxis')
title('Spektrogram Sound one')
subplot(2,3,6)
spectrogram(sound_two,NFFT,0,NFFT,fs,'yaxis')
title('Spektrogram Sound two')

    function callback_play(hObject,callbackdata,sig,fs)
       if max(abs(sig)) >= 1
             warning('Clipping detected. Using soundsc');
             soundsc(sig,fs);
       else
            soundsc(sig,fs);
       end
    end

end

