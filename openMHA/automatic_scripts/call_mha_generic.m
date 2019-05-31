% call_mha_generic(mha_config_file,infile,outfile,mha_dir,varargin)
%
% This function is used for calling the MHA using the configuration saved 
% in mha_config_file. It will process the wav-file infile and write the
% processed output to the wav-file outfile. The MHA starting script is located in
% mha_dir.
% 
% Optionally, you can specify different parameters, to e.g. change the configuration
% of a MHA plugin with a user-specific setting, like often done for hearing aid gains.
% The available parameters are described in the input parser in the code. 
% Example:
%
% mha_config_file    = 'my_mha.cfg';
% infile = 'MHA-tmp-in.wav';
% outfile = 'MHA-tmp-out.wav';
% % write signal to disk:
% audiowrite(infile, signal, fs, 'BitsPerSample', 32); % 32 Bits per Sample is important!
% % run MHA processing
% call_mha(mha_config_file,infile,outfile);
% % read in resulting signal
% [mha_sig,fs] = audioread(outfile);
% Example for modifying a mha plugin:
% call_mha(mha_config_file,infile,outfile,'plugin_config_file','myGainTable','plugin_name', 'mha.overlapadd.mhachain.dc')
% Authors: Alina Ernst
%          Ayham Zedan <ayham.zedan@uni-oldenburg.de>
%          Ben Williges <ben.williges@uni-oldenburg.de>

function call_mha_generic(mha_config_file,infile,outfile,varargin)

% Parse argument inputs
p = inputParser;
p.addRequired('mha_config_file',@ischar);
p.addRequired('infile',@ischar);
p.addRequired('outfile',@ischar);
p.addParameter('mha_dir','',@ischar); % Configuration file for the mha plugin you want to change
p.addParameter('plugin_config_file','',@ischar); % Configuration file for the mha plugin you want to change
p.addParameter('plugin_name','',@ischar); %String containing the MHA-specific full path to the plugin you want to change, e.g. 'mha.overlapadd.mhachain.dc'
p.addParameter('host','localhost',@ischar); %host name of the pc, where the MHA is started
p.addParameter('port',33337,@isnumeric); %Port, where the mha listen
% Validation
p.parse(mha_config_file, infile, outfile, varargin{:});
par = p.Results;

% ensure openmha matlab files are accessible
res = which('mha_query');
if strcmp(res,'')
    error('Could not find openmhas Matlab scripts in MATLAB PATH. Please add the directory mha/tools/mfiles to MATLAB PATH');
end
%start mha

if isunix
    if isempty(strfind(getenv('LD_LIBRARY_PATH'),'openMHA'))
        setenv('LD_LIBRARY_PATH', [':' par.mha_dir]);
    end
    
    if isempty(strfind(getenv('PATH'),'openMHA'))
        setenv('PATH', [getenv('PATH') ':' par.mha_dir]);
    end
%     [~,PID] = system(['xterm -hold -e ' 'mha --log=CON --port=' num2str(par.port) ' & echo $! ']);
%    
    [~,PID] = system(['mha --log=CON --port=' num2str(par.port) '&' ]);
else %For windows
system('run_MHA.cmd&')
end

pause(1)
work.mha = struct('host',par.host,'port',par.port);

mha_query(work.mha,'',['read:' par.mha_config_file]); 
mha_set(work.mha, 'io.in', infile);
mha_set(work.mha, 'io.out', outfile);
    
% Basic calibration:
mha_set(work.mha,' mha.calib_in.peaklevel', 120.0);
mha_set(work.mha, 'mha.calib_out.peaklevel', 120.0);
    
% Prepare MHA to run:
mha_set(work.mha,'cmd','prepare'); %MHA must do this at least once, otherwise the mhaplugins will have not the correct settings/default settings

% Correct plugin setup:
if ~isempty(par.plugin_config_file)
        mha_query(work.mha, par.plugin_name, ['read:' par.plugin_config_file]);
        mha_set(work.mha,'cmd','prepare'); %Now MHA should be ready to run
end    
% Run mha processing
mha_set(work.mha, 'cmd', 'start'); %Process files
mha_set(work.mha, 'cmd', 'quit'); %instead of release; MHA must be closed
if isunix
    kill_str = ['kill ' num2str(PID)]; 
    else %For windows
    kill_str = '"C:\Windows\System32\taskkill.exe" /F /im cmd.exe &';
end
[status1]=system('exit'); %This should close the cmd windows
pause(1)
end 

