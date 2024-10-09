%% Set paths for reconstruction

% check if luna or flux environment, set pathnames
hostname = getenv('HOSTNAME');
baseDir = '/opt/aumc-apps/';
if strcmp(hostname,'luna-01')
    bartpath = [baseDir 'bart/0.8.00-gpu/'];
    mreconpath = [baseDir 'matlab/toolbox/MRecon-5.4.2-mod/'];
    addpath('/opt/aumc-apps/matlab/toolbox/nifti_20140122');
else    % flux
    bartpath = [baseDir 'bart/0.7.00-gpu/'];
    mreconpath = [baseDir 'matlab/toolbox/MRecon-5.4.2-mod/'];
    addpath('/opt/amc/matlab/toolbox/nifti_20140122/');
end

%% BART SETTINGS
% addpath(fullfile('/opt/amc/bart-0.4.03', 'matlab'));
% setenv('TOOLBOX_PATH', '/opt/amc/bart-0.4.03/bin');
% /opt/aumc-apps/bart/0.8.00-gpu/ % luna, since 2023-12-19
% /opt/aumc-apps/bart/0.7.00-gpu/ % luna, since 2023-12-19

addpath(genpath(bartpath))
setenv('TOOLBOX_PATH', fullfile(bartpath,'bin/'));
setenv('OMP_NUM_THREADS','4');  % can be overwritten in .pro.dat

% scanner release / MRecon compatability
% ----scanner--------|--MRecon (from)------(up to)----
% 3T R5.1.8 SWID31 :    ?               ?
% 3T R5.3.1 SWID57 :    MRecon-3.0.541  MRecon-3.0.541
% 3T R5.4   SWID420:    MRecon-3.0.545
% 7T R5-B   SWID129:    MRecon-3.0.529  MRecon-3.0.537

% mreconpath = '/opt/aumc-apps/matlab/toolbox/MRecon-5.4.1/';   % since 2023-12-19
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
% mreconpath = '/opt/amc/matlab/toolbox/MRecon-4.3.1-mod/'; % ~ since 2022-09-08 for Release R7 --> does not work for release R9 or higher
% mreconpath = [baseDir 'matlab/toolbox/MRecon-5.4.1/']; % since 2023-12-21
% mreconpath = [baseDir 'matlab/toolbox/MRecon-5.4.2-mod/']; % since 2024-03-22
addpath(genpath(fullfile(mreconpath)))

%% enable warnings
warning('on','all')
