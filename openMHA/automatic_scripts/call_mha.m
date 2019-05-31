function call_mha(mha_config_file,infile,outfile)
portNum = 33337;
% [xtermstatus,~] =  system('pgrep -x "xterm"');
% if xtermstatus
    [~,PID] = system(['xterm -hold -e ' '../' filesep 'openMHA' filesep 'automatic_scrtipt_MHA_Bash ' num2str(portNum) ' & echo $! ']);
     pause (1.5)
% end
work.mha = struct('host','localhost','port',portNum);
mha_set(work.mha, 'io.in', infile);
mha_set(work.mha, 'io.out', outfile);

% basic calibration:
mha_set(work.mha,' mha.calib_in.peaklevel', 120.0);
mha_set(work.mha, 'mha.calib_out.peaklevel', 120.0);

% prepare MHA to run:
mha_set(work.mha,'cmd','prepare'); %MHA must do this at least once, otherwise the mhaplugins will have not the correct settings/default settings
%         mha_query(work.mha,'','save:log.cfg');

% correct fitting:
mha_query(work.mha, 'mha.overlapadd.mhachain.dc', ['read:' mha_config_file]); %This file contains the hearing loss depending gains.
% Obtain the correct file by:
% 1) Start MHA, e.g. by running extern_MHA('MHA_output_to_wav')
% with a breakpoint after the 'prepare' cmd has been issued (needed
% for correct channel allocations)
% 2) (in Matlab-subdirectory): mhagui_fitting
% 3) Enter Audiogramm, then Click ok
% 4) create first fit
% 5) in MHA fitting tool, choose camfit_cicand_left
% 6) Close mhacontrol-window, run mhagui_generic
% 7) Change to configuration of dynamic compressor (here:
% mha.mhachain.ola.sg.mbc.dc and save this configuration under the
% name you later want to set via def.mha_config_file
% Alternative: mha_query(work.mha, 'mha.mhachain.ola.sg.mbc.dc', 'save:Example_fitting_NALNL1.cfg')


%     mha_query(work.mha, 'mha.overlapadd.mhachain.dc', 'save:Example_fitting_NALNL1.cfg')
mha_set(work.mha,'cmd','prepare'); %Now MHA should be ready to run

% run mha processing
mha_set(work.mha, 'cmd', 'start'); %Process files
mha_set(work.mha, 'cmd', 'release');
mha_query(work.mha,'','save:log_config_file.cfg');
%     mha_set(work.mha, 'cmd', 'quit');

system (['kill ' num2str(PID)]);
end