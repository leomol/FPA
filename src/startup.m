% Add path to dependencies.
setup()
function setup()
    root = fileparts(mfilename('fullpath'));
    addpath(root);
    addpath(genpath(root));
end