function [Formants]= FormantEstimator(Fs, vFrame)
% retreived form http://de.mathworks.com/help/signal/ug/formant-estimation-with-lpc-coefficients.html
Formants=zeros(1,12);
n=18; % Order for LPC
LPCoeff = lpc(vFrame,n);
NoFormants=find(isnan(LPCoeff));
if length(NoFormants) < 1
    rts = roots(LPCoeff );
    rts = rts(imag(rts)>=0);
    angz = atan2(imag(rts),real(rts));
    [frqs,indices] = sort(angz.*(Fs/(2*pi)));
    bw = -1/2*(Fs/(2*pi))*log(abs(rts(indices)));
    nn = 1;
    VarifiedFormants=find(bw <400);
    NVarifiedFormants=length(find(frqs(VarifiedFormants)>150)); % formant freq should be larger than 150 Hz
    for kk = 1:NVarifiedFormants
        if (frqs(kk) > 150 && bw(kk) <400) % doppelt gemoppelt
            Formants(nn) = frqs(kk);
            if nn == 1 && Formants(nn) < 900 && Formants(nn) > 150 % limit range of F1
                nn = nn+1;
            else if nn == 2 && Formants(nn) < 2500 && Formants(nn) > 450 % limit range of F2
                    nn = nn+1;
                else
                    Formants(nn)=[];
                end
           end
        end
    end
end
Formants=Formants(1:2);


% Alternative
% http://iitg.vlab.co.in/?sub=59&brch=164&sim=615&cnt=2
% https://www.clear.rice.edu/elec431/projects96/digitalbb/fmnts.html
% https://ccrma.stanford.edu/~jos/fp/Formant_Filtering_Example.html