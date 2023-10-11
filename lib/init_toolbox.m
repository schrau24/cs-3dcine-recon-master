%% Set MATLAB path and TOOLBOX_PATH environment variable (BART)
% This script is equal to the '/opt/amc/bart/vars.m' function


%% BART SETTINGS, insert bart path
% addpath(fullfile('/PATH/TO/BART/','matlab'))
% setenv('TOOLBOX_PATH', '/PATH/TO/BART/bin');
addpath(fullfile('/opt/amc/bart-0.5.00-gpu', 'matlab'));
setenv('TOOLBOX_PATH', '/opt/amc/bart-0.5.00-gpu/bin');

% for parallel computing
setenv('OMP_NUM_THREADS','4');

%% MRECON SETTINGS
% mreconpath = '/PATH/TO/MRECON';
mreconpath = '/opt/amc/matlab/toolbox/MRecon-4.3.1-mod/'; % default (latest)

addpath(genpath(fullfile(mreconpath)))

%% enable warnings
warning('on','all')