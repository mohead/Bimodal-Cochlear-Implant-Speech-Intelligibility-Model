function [parameters] = freq_distortion_params(varargin)

% Author: R. Bennett, 2011, modified by B. Williges, 2016

% ## PARAMETERWAHL
% Input parser (make all parameters optional)
parameters = inputParser; %creates Input parser structure
% ## Anzahl der Iterationen des Algorithmus
parameters.addParamValue('numIterations', 10,@isnumeric); % Minimum 2

% ## Bandwidthfactor (Bandbreitenfaktor)
parameters.addParamValue('bandwidthFactor', 4,@isnumeric); %  Default 4

% ## Bandwidthfactor only for Normalhearing Simulation
% ## NH-Bandbreitenfaktor
% ## Werte größer 1 für Normalhörendensimulation 
% ## Bandbreitenfaktor (oben) muss dann gleich 1 gesetzt werden.
parameters.addParamValue('bandwidthFactor_NH', 1,@isnumeric);
%--------------------------------------
% ## Samplingfrequenz
parameters.addParamValue('sampleFreqHz', 44100,@isnumeric);

% ## HI Filterordnung Algorithmus und Modell
parameters.addParamValue('gammaOrderHI', 4,@isnumeric);

% ## NH Gesamtfiterodnung
parameters.addParamValue('gammaOrderNH12', 4,@isnumeric); % Default 4

% ## NH-Ordnung 1
parameters.addParamValue('splitOrder1', 3,@isnumeric); % Default 3
%%## NH-Ordnung 2
parameters.addParamValue('splitOrder2', 1,@isnumeric); % Default 1 (4-3 = 1)

% ## Order of final reference
% ## NH Ordnung Modell 
parameters.addParamValue('gammaOrderNH', 4,@isnumeric);

% ## Filter density on ERB_aud Scale
% ## ERB-Faktor
parameters.addParamValue('filtersPerERBaud', 1,@isnumeric);

% ## Filterbankkonfiguration
% ## Basefrequency
parameters.addParamValue('baseFreqHz', 500,@isnumeric); % 500 Hz

% ## Filterbankkonfiguration
% ## Untere Grenzfrequenz
parameters.addParamValue('lowerCutoffFreqHz', 80,@isnumeric);	 % Achtung: Wenn Fehlermeldung, dann diesen Wert höher setzen.

% ## Obere Grenzfrequenz
parameters.addParamValue('upperCutoffFreqHz', 10000,@isnumeric)% 10000 Achtung: Wenn Fehlermeldung, dann diesen Wert tiefer setzen.

% ## Filter Delays
parameters.addParamValue('filterDelaySec', [0.016 0.016],@isnumeric);
parameters.addParamValue('debug',0,@isnumeric); %can be 0 (= no debug), 1 (= display text information) or 2(=additionally plot outputs)
% Validation
parameters.parse(varargin{:});
parameters = parameters.Results;
% ## NH-Ordnung 2
assert(isequal(parameters.splitOrder2,parameters.gammaOrderNH12 - parameters.splitOrder1),'splitOrder2 is not set correctly');
