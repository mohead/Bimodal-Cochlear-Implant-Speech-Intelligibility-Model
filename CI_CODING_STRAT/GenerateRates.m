function  ExSig= GenerateRates(par,mFormantData,FrameS)
% Function to generate a template of zeros and ones that is multiplied with
% pulse cycles later. Ones occur with respect to F0, F1 and 900 pps
%
% Input:
% --------
%         par             = Parameters set in an external parameter file
%         mFormantData    = [Time; VoicedUnvoiced; F0, F1; F2] interpolated
%         FrameS          = Frame Shift
% Output:
% --------
%         ExSig           = Rate Templates (Matrix)
%         
%
% Anja Eichenauer, 05.02.2016
% Universität Oldenburg


BaseRate=900; % Rate for F2 channels
mFormantData(5,:)=BaseRate; % Former F2 Frequency is changed to base rate
nFrames=size(mFormantData,2);
% generate Periodicity
vTPeriod=floor(1./(mFormantData(3:5,:)+eps).*((par.CIFs/FrameS)*par.n)); % every x Samples 1s are generated else 0, devided by FrameS due to block processing, times par.n because during one period n electrode slots have to be cycled
vTPeriodUnvoiced=(1/BaseRate)*((par.CIFs/FrameS)*par.n);
ExSig=zeros(size(vTPeriod,1),nFrames); % preallocation
LastIDX=ones(size(vTPeriod,1),1); % last idx of previous period
ExponentialDecay(1:par.n)=1; % used to be exponential decay, now single pulses
ExponentialDecayUnvoiced(1:par.n)=1;

for uu=1:nFrames % over time
    for zz=1:size(vTPeriod,1) % over Formants
        if mFormantData(2,uu) == 1 && vTPeriod(zz,uu)<= floor(1./80.*((par.CIFs/FrameS)*par.n)) && vTPeriod(zz,uu)>= floor(1./900.*((par.CIFs/FrameS)*par.n)) % Muss als voiced eingestuft sein und f zw. 80 Hz u. 900 Hz
            
            if uu >= LastIDX(zz)% is loop larger than previous idx of period
                PeriodIDX=LastIDX(zz)+1:LastIDX(zz)+vTPeriod(zz,uu); % all idxs of current period
                ExSig(zz,PeriodIDX(1):PeriodIDX(1)+size(ExponentialDecay,2)-1)=ExponentialDecay; % 12 ones for one whole cycle, later the matching channels are allocated
                LastIDX(zz)=PeriodIDX(end); % update last idx of current period
            end
            
        else % unvoiced case
            if uu >= LastIDX(zz)
                PeriodIDX=LastIDX(zz)+1:LastIDX(zz)+vTPeriodUnvoiced; % all idxs of current period
                ExSig(zz,PeriodIDX(1):PeriodIDX(1)+size(ExponentialDecayUnvoiced,2)-1)=ExponentialDecayUnvoiced; % 12 ones for one whole cycle
                LastIDX(zz)=PeriodIDX(end); % update last idx of current period
            end
        end
    end
end