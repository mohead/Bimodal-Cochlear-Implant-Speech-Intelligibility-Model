function cellArgs	= nargdef(FuncVarargin, varargin)
%
% function cellArgs	= nargdef(FuncVarargin, [ListOfDefArguments])
%
% This function replaces missing or empty input 
% arguments to a function. 
%
% Assuming a function declared like
%    x = myfunc(varargin)
% which should have 3 input Arguments like in
%    x = myfunc(X,Y,Z)
% may has default input values like
%    X = [1:10]
%    Y = []
%    Z = 'string'
% then the following sequence patches the default 
% values which are empty or missing in varargin:
%    cellArgs	= nargdef(varargin, [1:10], [], 'string')
%    X = cellArgs{1};
%    Y = cellArgs{2};
%    Z = cellArgs{3};
%

% fill defaults
cellArgs	= varargin;

% replace with given arguments
for i = 1:length(FuncVarargin),
	if ~isempty(FuncVarargin{i})
		cellArgs{i}	= FuncVarargin{i};
	end;
end;

return;
