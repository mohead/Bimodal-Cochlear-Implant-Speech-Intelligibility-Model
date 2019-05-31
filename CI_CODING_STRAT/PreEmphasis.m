function vPreEmp = PreEmphasis(Input)

%%-- Stefan Fredelakes Preemp
fc    = 1200;
w     = 2*fc/16000; % fc/fs/2 = 2*fc/fs = 1200 Hz
[b,a] = butter(1,w,'high');
vPreEmp = filter(b,a,Input);
