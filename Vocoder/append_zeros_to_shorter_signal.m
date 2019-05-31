function [sig1,sig2] = append_zeros_to_shorter_signal(sig1,sig2)
% This function compares the length of the two input signals and will
% append zeros to the shorter one, so that the two signals will match in
% length
assert(isequal(length(sig1),size(sig1,1)),'Sig1 must be of dimension samples x channels');
assert(isequal(length(sig2),size(sig2,1)),'Sig2 must be of dimension samples x channels');

if length(sig1)>length(sig2)
    add_zeros = zeros(length(sig1)-length(sig2),size(sig2,2));
    sig2 = [sig2; add_zeros];
elseif length(sig1)<length(sig2)
    add_zeros = zeros(length(sig2)-length(sig1),size(sig1,2));
    sig1 = [sig1; add_zeros];
else 
    %same length for both sig
end
assert(isequal(length(sig1),length(sig2)),'Error when adjusting signal lengths');
end