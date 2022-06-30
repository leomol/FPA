function setup()
    root = fileparts(mfilename('fullpath'));
    addpath(root);
    addpath(fullfile(root, 'common'));
    addpath(fullfile(root, 'common/TDTbin2mat'));
end