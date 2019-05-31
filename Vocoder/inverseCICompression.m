% inverse function to CI compression

function out = inverseCICompression(in,B,M,alpha_c)
    % loudness inversion:
    out = (exp(log(1+alpha_c).*in)-1)./alpha_c .*(M-B)+B;
    out(in==0) = 0; %remove DC-offset introduced by adding B