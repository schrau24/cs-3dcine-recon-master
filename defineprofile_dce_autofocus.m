%% define recon profile
% Edit and run this script to define a recon profile

% (!) Disabled parameter (commented) will be filled in automatically

% (!) Mandatory parameters need to be filled in!
%% housekeeping / initialization
restoredefaultpath;
addpath(fullfile(cd,'./lib/'));
% addpath(genpath(fullfile(cd,'./Elastix/')));
clear;
close all;
clc;

%% MANDATORY PARAMETERS
profile = {};

% have the user grab the data_dir and data_target, the coilsurvery and
% senseref are automatically determined
homeDir = fileparts(mfilename('fullpath'));

% define folders and files
baseDir = ['/home/' getenv('USER') '/lood_storage/divi/Projects/asap/'];
cd(baseDir)
[filename, pathname] = uigetfile('*.raw','Select raw data file');
cd(homeDir)
% use the folder name as the base for the profile name
s = filename(1:end-4);
% disp(['profile name: ' s])
profile.name            = [s '_DCE_autofocus'];

profile.data_dir        = pathname;
profile.data_target     = filename(1:end-4);
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
profile.bart_pics_cmd       = 'pics -R T:2048:0:0.01 -R T:7:0:0.001 -i 40 -S -d5';
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
% profile.Recon_SENSE                         = 'Yes';
% profile.Recon_DcOffsetCorrection            = 'Yes';
% profile.Recon_PDACorrection                 = 'Yes';
% profile.Recon_RandomPhaseCorrection         = 'Yes';
% profile.Recon_MeasPhaseCorrection           = 'Yes';
% profile.Recon_PartialFourier                = 'Yes';
% profile.Recon_Gridding                      = 'Yes';
% profile.Recon_RingingFilter                 = 'Yes';
% profile.Recon_RingingFilterStrength         = [0.25, 0.25, 0.25];
% profile.Recon_EPIPhaseCorrection            = 'Yes';
% profile.Recon_EPICorrectionMethod           = 'Linear';
% profile.Recon_EPI2DCorr                     = 'Yes';
% profile.Recon_EPICorrPerLocation            = 'No';
% profile.Recon_RotateImage                   = 'Yes';
% profile.Recon_GeometryCorrection            = 'Yes';
% profile.Recon_RemoveMOversampling           = 'Yes';
% profile.Recon_RemovePOversampling           = 'Yes';
% profile.Recon_ConcomitantFieldCorrection    = 'Yes';
% profile.Recon_FlowPhaseCorrection           = 'Yes';
% profile.Recon_DivideFlowSegments            = 'Yes';
% profile.Recon_TKE                           = 'No';
% profile.Recon_Venc                          = [150, 150, 150];
% profile.Recon_kv                            = ''; % is array
% profile.Recon_FluidDensity                  = 1060;
% profile.Recon_Average                       = 'Yes';
% profile.Recon_ArrayCompression              = 'No';
% profile.Recon_ACNrVirtualChannels           = '';
% profile.Recon_ACMatrix                      = '';
% profile.Recon_ImmediateAveraging            = 'Yes';
% profile.Recon_ExportRECImgTypes             = 'M';
% profile.Recon_AutoUpdateInfoPars            = 'Yes';
% profile.Recon_Sensitivities                 = '';
% profile.Recon_SENSEPsi                      = '';
% profile.Recon_SENSERegStrength              = 2;
% profile.Recon_StatusMessage                 = 'Yes';
% profile.Recon_Logging                       = 'No';
% profile.Recon_AutoChunkHandling             = 'Yes';
% profile.Recon_EddyCurrentCorrection         = 'No';

% profile.reco_skip_PW                        = 0; % skip Pre-whitening        [ 1:yes, 0:no]
% profile.reco_skip_CLEAR                     = 0; % skip CLEAR correction     [ 1:yes, 0:no]
profile.reco_PCMRA                          = 0; % Create PC-MRA             [ 1:yes, 0:no]
profile.reco_NoiseClip  = 'No';
profile.reco_NoiseClipValue     = 100;

%% % ENCODING PARAMETERS
% profile.Encoding_KzOversampling             = 1; %

%% % DCE PARAMETERS
profile.DCETimeResolution = 3;              % in seconds

%% % RESPIRATORY CORRECTION PARAMS
profile.respCorrection = 1;
profile.respPhases = 8;
%% % Retrospective Undersampling
% profile.RU_do   = 0; % 1: do retrospective undersampling
% profile.RU_mask = ''; % file location needed, if RU_do = 1

%% save profile
disp(['DCE time resolution = ' num2str(profile.DCETimeResolution) ' seconds'])
disp('Saving reconstruction profile:')
disp([profile.name '.pro.dat'])
save_profile(profile);
