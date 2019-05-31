function [mFormantData,  Output, FormantDataInPl,NewOutput, NewTime, StimOrder, vDistance,AllFormantCh] = FormantLocking(vSignal, par, vTimeAll,OuterCornerFrqs,FL_ChannelAllocation, EnvSig, FrameS, nSlots)


[mFormantData, FormantDataInPl] = ExtractFormantData(vSignal, par,vTimeAll);

ExSig = GeneratePeriodicity(par,mFormantData,FrameS); % Generate Modulation Shape (ExSig = Excitation signal for signal modulation)
[Output, NewOutput, NewTime, StimOrder, vDistance,AllFormantCh] = AllocateFormantMod2Ch(OuterCornerFrqs, FL_ChannelAllocation,EnvSig,ExSig,mFormantData,par, nSlots);

