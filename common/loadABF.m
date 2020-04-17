% [data, units, names] = loadABF(filename)
% Load abf files from Axon.
% 
% This is currently a function draft and may only work for special cases
% (e.g. one channel and gap-free mode).

% 2020-03-02. Leonardo Molina.
% 2020-03-09. Last modified.
function [data, units, names, header] = loadABF(filename)
    [data, ~, header] = resources.abfload(filename, 'doDispInfo', false);
    units = header.recChUnits;
    names = ['time'; header.recChNames];
    
    duration = header.lActualAcqLength * header.fADCSampleInterval * 1e-6;
    nSamples = size(data, 1);
    time = linspace(0, duration, nSamples);
    data = [time(:), data];
end