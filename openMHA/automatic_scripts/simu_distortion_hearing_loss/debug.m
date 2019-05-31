%% function debug(msg,flags,fun_handle,varargin)
% debug function to help debug matlab programms/scripts.
% Depending on the state of flags.debug this function will output nothing
% (flags.debug = 0), only display a message-string 'msg' (flags.debug = 1),
% or execute the function handle 'fun_handle' with additional parameters
% provided by varargin. This function might be useful for debugging
% algorithms or models. Just specify a debug-flag inside of your
% parameter/flag structure and then write a line of code for each stage, where
% you might want to display additional informations. Thus it is very easy
% to switch between debug-mode and "running" stage of your programm,
% without having to comment out huge amounts of codes or placing if/end all
% over the code.
% 
% Example usage:
% aa = [1:10];
% flags.debug = 1; %This will just show the message-part
% bb = aa+2;
% debug('Finished processing of bb',flags);
% flags.debug = 2; %This will also evaluate the function handle
% debug('Finished processing of bb',flags,@plot,bb,'r'); %Additionally plot
% the output
function debug(msg,flags,fun_handle,varargin)
if nargin == 2 % If you don't provide more than two inputs, just display the message
    flags.debug = 1;
end
if flags.debug == 1 % Just show simple debug-message
    disp(msg); % Prints msg to screen.
elseif flags.debug == 2
    disp(msg); %Prints msg to screen
    fun_handle(varargin{:}) %Evaluates function handle with additional input using varargin
end
end