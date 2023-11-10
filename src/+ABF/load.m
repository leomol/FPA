% [data, units, names] = ABF.load(filename)
% Load abf files from Axon.
% 
% This is currently a function draft and may only work for special cases
% (e.g. one channel and gap-free mode).

% 2020-03-02. Leonardo Molina.
% 2023-11-09. Last modified.
function [data, units, names, header] = load(filename)
    [data, ~, header] = abfload(filename, 'doDispInfo', false);
    units = header.recChUnits;
    names = ['time'; header.recChNames];
    
    duration = header.lActualAcqLength * header.fADCSampleInterval * 1e-6;
    nSamples = size(data, 1);
    time = linspace(0, duration, nSamples);
    data = [time(:), data];
end