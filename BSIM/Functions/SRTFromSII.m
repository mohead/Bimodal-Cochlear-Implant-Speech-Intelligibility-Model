function [SRT,slope] = SRTFromSII(E,N,T,I,G,fRefLev,a,b)

RefSII = b;

L    = -64;
fSII =   0;
vSteps = 2.^[4:-2:-6];
iSign  = +1;

for step = 1:length(vSteps)
    while fSII*iSign < RefSII*iSign
        fSII = SII(E+L,N,T,I,G);
        L = L + iSign * vSteps(step);
    end
    iSign = iSign * -1;
end

SRT = L;
