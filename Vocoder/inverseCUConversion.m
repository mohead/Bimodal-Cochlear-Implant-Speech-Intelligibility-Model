function out = inverseCUConversion(in,TCL, MCL, Volume)
    % Inverse Conversion from CU to acoustic:
    acoustic = bsxfun(@minus, in, TCL);
    for ii = 1:size(acoustic,1)
         idx = (acoustic(ii,:) == -TCL(ii));
         acoustic(ii,idx) = 0; % Prevent overshoot
    end
    out = bsxfun(@rdivide, acoustic, Volume*(MCL-TCL));