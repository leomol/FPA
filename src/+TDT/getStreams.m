% streams = getStreams(folder)
% Parses a project folder recorded with TDT DAQ / Synapse and returns the
% name of the available streams.

% 2019-02-01. Leonardo Molina.
% 2025-11-24. Last modified.
function streams = getStreams(folder)
    data = TDTbin2mat(folder, 'T1', 0, 'T2', 1);
    streams = fieldnames(data.streams)';
end