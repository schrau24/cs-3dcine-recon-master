%% define recon profile
% Edit and run this script to define a recon profile

% (!) Disabled parameter (commented) will be filled in automatically

% (!) Mandatory parameters need to be filled in!
%% housekeeping / initialization
restoredefaultpath;
addpath(fullfile(cd,'./lib/'));
clear;
close all;
clc;

%% MANDATORY PARAMETERS
profile = {};

% have the user grab the data_dir and data_target, the coilsurvery and
% senseref are automatically determined

% define folders and files
% baseDir = ['/home/' getenv('USER') '/lood_storage/divi/Projects/dcerecon/Spirals_for_DCE/'];
baseDir = ['/home/' getenv('USER') '/lood_storage/divi/Projects/asap/WP3/Data/WP3b_patients/'];

[filename, pathname] =  uigetfile([baseDir '*.raw'],'Select one of the VFA raw data files');

% now find all the vfa files, note 'vfa' needs to be in the name
allFiles = dir([pathname '/*vfa*.raw']);

% now just loop over them and make profiles
for ff = 1:length(allFiles)
    
    % use the folder name as the base for the profile name
    tmp = allFiles(ff).name;
    s = tmp(1:end-30);
    % disp(['profile name: ' s])
    profile.name            = [s '_t1map_autofocus'];
    
    profile.data_dir        = pathname;
    profile.data_target     = tmp(1:end-4);
    s = dir([pathname '*senseref*']);
    if ~isempty(s)
        s = s(1).name;
        profile.data_senseref   = s(1:end-4);
    else
        profile.data_senseref = 'noSenseRefScan';
        disp('No senseref scan found, returning empty')
    end
    
    s = dir([pathname '*coilsurvey*']);
    if ~isempty(s)
        s = s(1).name;
        profile.data_coilsurvey = s(1:end-4);
    else
        profile.data_coilsurvey = 'noCoilSurveyScan';
        disp('No coilsurveyscan found, returning empty')
    end
    
    %% % OPTIONAL PARAMETERS
    
    %% define export [ 1:yes, 0:no]
    profile.exp_mask        = 0;
    profile.exp_parrec      = 1;
    profile.exp_nii         = 1;
    profile.exp_dcm         = 0;
    profile.exp_sense       = 0;
    
    %% % BART: Compressed Sensing
    %  % [ 1:yes, 0:no]
    profile.bart_pics           = 1;
    % 1: use BART PICS ADMM reconstruction
    % 0: use MRecon reconstruction
    profile.bart_sensemap       = 1;
    % 1: use estimated bart sense maps
    % 0: use measured sense maps (MRSense)
    profile.bart_cc             = 0;
    % 1: use bart coil compression
    % 0: no coil compression
    profile.bart_skip_fmac_SENSEUnfold = 1;
    % 1: skip fmac after pics and also skip SENSEUnfold in MRecon
    % 0: use fmac after pics and also SENSEUnfold in MRecon
    profile.bart_pics_cmd       = 'pics -R T:7:0:0.01 -i 40 -S -d5';
    % pics [-l ...] [-r f] [-c] [-s f] [-i d] [-t <string>] [-n] [-g] [-p <string>] [-I ...] [-b d] [-e] [-W <string>] [-d d] [-u f] [-C d] [-f f] [-m ...] [-w f] [-S] [-B d] [-K] <kspace> <sensitivities> <output>
    
    profile.bart_sensemap_cmd   = 'caldir 20';
    % profile.bart_sensemap_cmd   = 'ecalib -m1 -I';
    % ecalib [-t f] [-c f] [-k ...] [-r ...] [-m d] [-S] [-W] [-I] [-1] [-v f] [-a]
    % caldir cal_size
    
    profile.bart_cc_cmd         = 'cc -p 8 -E';
    % cc [-p d] [-M] [-r ...] [-A] [-S ...] [-G ...] [-E ...] <kspace> <coeff>|<proj_kspace>
    % Types - S: SVD; G: Geometric ;E: ESPIRiT
    
    profile.bart_parpool             = 12;
    % 0:     NO  parallel computing in CS reconstruction
    % N > 0: Use parallel computing in CS reconstruction with N workers
    
    %% % RECON PARAMETER
    profile.Recon_CoilCombination               = 'pc';
    profile.Recon_ImageSpaceZeroFill            = 'No';
    profile.Recon_kSpaceZeroFill                = 'No';
    
    profile.reco_PCMRA                          = 0; % Create PC-MRA             [ 1:yes, 0:no]
    profile.reco_NoiseClip  = 'No';
    profile.reco_NoiseClipValue     = 100;
    
    %% % DCE PARAMETERS
    profile.DCETimeResolution = 60;              % in seconds
    
    %% % RESPIRATORY CORRECTION PARAMS
    profile.respCorrection = 1;
    profile.respPhases = 8;

    %% save profile
    disp(['DCE time resolution = ' num2str(profile.DCETimeResolution) ' seconds'])
    disp('Saving reconstruction profile:')
    disp([profile.name '.pro.dat'])
    save_profile(profile);
end