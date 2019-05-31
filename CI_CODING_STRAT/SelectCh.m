function nChSig = SelectCh(vSigCh,par)
%% function to select n out of m channels with highest amplitude
% Additionally, deactivated channels are sorted out here.
% Input:  
% --------
%         vSigCh    = Weighted Signal envelope in Channels
%         par       = Parameters set in an external parameter file
% Output: 
% --------
%         nChSig    = Weighted Signal envelope in Channels, 0-amplitude if
%                     sorted out
%         

nChSig=zeros(par.m,1);
vSigCh(par.nDeactivated)=-Inf; % set amplitudes of deactivated ch -Inf
vSigCh(par.pps==0)=-Inf; % if this should ever happen!?!?!
[~, Idx] = sort(vSigCh,'descend'); % sort channels with respect to amplitude
nChSig(Idx(1:par.n))=vSigCh(Idx(1:par.n)); % keep n highest amplitudes, others stay 0
