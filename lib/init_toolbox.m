%% Set MATLAB path and TOOLBOX_PATH environment variable (BART)
% This script is equal to the '/opt/amc/bart/vars.m' function


%% BART SETTINGS
% addpath(fullfile('/opt/amc/bart-0.4.03', 'matlab'));
% setenv('TOOLBOX_PATH', '/opt/amc/bart-0.4.03/bin');
addpath(fullfile('/opt/amc/bart-0.5.00-gpu', 'matlab'));
setenv('TOOLBOX_PATH', '/opt/amc/bart-0.5.00-gpu/bin');
setenv('OMP_NUM_THREADS','4');

%% MRECON SETTINGS
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-4.3.1-mod/'; % default (latest)
mreconpath = '/opt/amc/matlab/toolbox/MRecon-5.4.2/';		% latest version, 
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-5.1.0/';       	% 202309 to 20230426, for .sin file reading (fetal)
% mreconpath = '/opt/amc/matlab/toolbox/MRecon'; % default (latest)

% scanner release / MRecon compatability
% ----scanner--------|--MRecon (from)------(up to)----
% 3T R5.1.8 SWID31 :    ?               ?
% 3T R5.3.1 SWID57 :    MRecon-3.0.541  MRecon-3.0.541
% 3T R5.4   SWID420:    MRecon-3.0.545
% 7T R5-B   SWID129:    MRecon-3.0.529  MRecon-3.0.537

% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.482'; % since 2015-04-08
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.506'; % since 2015-10-21
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.515'; % since 2016-01-19
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.519'; % since 2016-03-22
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.523'; % since 2017-01-06
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.529'; % since 2017-01-06
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.532'; % since 2017-03-21
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.535'; % since 2017-05-17
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.537'; % since 2017-06-27
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.539'; % since 2017-07-25
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.541'; % since 2017-08-09
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.545'; % since 2017-09-07
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.553'; % since 2018-06-26 -- buggy in DICOM export 
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.554'; % since 2018-06-07
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.556'; % since 2018-07-10
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-3.0.557'; % since 2018-09-27

addpath(genpath(fullfile(mreconpath)))

%% enable warnings
warning('on','all')
