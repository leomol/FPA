% 2021-02-25. Leonardo Molina.
% 2021-02-25. Last modified.

% Export data as csv files.
% inputDataFile = ...;
% configuration = ...;
% results = FPA(...);

[folder, basename] = fileparts(inputDataFile);
exportCSV(folder, basename, results, configuration);

function exportCSV(folder, basename, results, configuration)
    % Save data for post-processing.
    nConditions = numel(configuration.conditionEpochs) / 2;
    triggeredWindow = size(results.windowIds, 2);
    halfWindow = (triggeredWindow - 1) / 2;
    windowTemplate = -halfWindow:halfWindow;
    triggeredWindowHeader = [repmat('n, ', 1, halfWindow), 'zero', repmat(', p', 1, halfWindow)];
    triggeredWindowFormat = ['%i', repmat(', %.4f', 1, triggeredWindow), '\n'];

    % File #1: time vs dff.
    % Rows represent increasing values of time with corresponding dff values.
    output = fullfile(folder, sprintf('%s - dff.csv', basename));
    fid = fopen(output, 'w');
    fprintf(fid, '# time, dff\n');
    fprintf(fid, '%.4f, %.4f\n', [results.time, results.dff]');
    fclose(fid);
    
    % File #2: AUC.
    output = fullfile(folder, sprintf('%s - AUC.csv', basename));
    fid = fopen(output, 'w');
    fprintf(fid, '# condition, area, duration\n');
    fprintf(fid, '%i, %.4f, %d\n', [(1:nConditions)', results.area, results.duration]');
    fclose(fid);
    
    % File #3: triggered windows with corresponding epoch label.
    % Order depends on epoch definitions. Overlapping is possible.
    % Rows represent a single peak:
    % First column is the condition label of the peak and is followed by the trace around each peak, with each peak at the center column (n / 2 + 1) labeled with c.
    output = fullfile(folder, sprintf('%s - peaks.csv', basename));
    fid = fopen(output, 'w');
    fprintf(fid, '# condition, %s\n', triggeredWindowHeader);
    fprintf(fid, ['#          ', sprintf(', %.2f', windowTemplate / results.frequency), '\n']);
    fprintf(fid, triggeredWindowFormat, [results.windowLabels, results.dff(results.windowIds)]');
    fclose(fid);
    
    % File #4: Average of the above.
    output = fullfile(folder, sprintf('%s - average peaks.csv', basename));
    fid = fopen(output, 'w');
    fprintf(fid,  '# condition, %s\n', triggeredWindowHeader);
    fprintf(fid, ['#          ', sprintf(', %.2f', windowTemplate / results.frequency), '\n']);
    fprintf(fid, triggeredWindowFormat, [unique(results.windowLabels, 'stable'), results.windowDff]');
    fclose(fid);
end