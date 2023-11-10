% CSV.load
% Load data from csv files.

% 2023-11-09. Leonardo Molina
% 2023-11-09. Last modified
function varargout = load(varargin)
    varargout = cell(1, nargout);
    [varargout{:}] = loadData(varargin{:});
end