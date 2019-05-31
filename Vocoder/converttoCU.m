    function electrodogramm_CU = converttoCU(electrodogramm, TCL, MCL, Volume)
    % function electrodogramm_CU = converttoCU(electrodogramm, TCL, MCL, Volume)
    % 
    % This function transforms an electrodogramm amplitude to a fixed range of
    % Threshold (TCL) and maximum (MCL) Clinical Units [CU]. 
    % The mapping function used is described in page 68 of
    % Brett Swansons Dissertation "Pitch percept with cochlear implants",
    % equation 4.26 in
            assert(all(TCL < MCL), 'TCL must be smaller than MCL!');
            assert(all(MCL <= 1200), 'MCL must be smaller than 1200!');
            assert(Volume <=1 && Volume >=0, 'Volume must be between [0 1]');
            electrodogramm_CU = bsxfun(@times,(MCL-TCL), electrodogramm.*Volume);
            electrodogramm_CU = bsxfun(@plus, TCL, electrodogramm_CU);
            % Remove TCL-offset (i.e. silence (zeros) would be transformed
            % to the TCL-Current. This means that at all times the
            % CI-Patient would hear a tone, which is not realistic at all)
            for ii = 1:size(electrodogramm,1)
               idx = (electrodogramm_CU(ii,:) == TCL(ii));
               electrodogramm_CU(ii,idx) = 0; 
            end
    end