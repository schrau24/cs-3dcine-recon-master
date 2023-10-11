%% define recon profile
% Edit and run this script to define a recon profile

% (!) Disabled parameter (commented) will be filled in automatically

% (!) Mandatory parameters need to be filled in!

%% MANDATORY PARAMETERS
profile = {};

% define folders and files
profile.name            = 'name';
profile.data_dir        = ['/home/' getenv('USER') '/lood_storage/divi/Projects/PROJECTNAME/incoming/parrec/...'];
profile.data_target     = '....raw';
%profile.data_senseref   = 'vo_13102022_1829231_1000_19_wip_senserefscanV4.raw';
%profile.data_coilsurvey = 'vo_13102022_1829128_1000_16_wip_coilsurveyscanV4.raw';
profile.output_dir = ['/home/' getenv('USER') '/lood_storage/divi/Projects/PROJECTNAME/incoming/recons/SUBJECTNAME/3d_cine/'];

%% % OPTIONAL PARAMETERS

%% define export [ 1:yes, 0:no]
profile.exp_mask        = 1;
profile.exp_parrec      = 1;
profile.exp_dcm         = 1;
profile.exp_sense       = 1;
profile.exp_meanResp    = 1;
profile.exp_mat         = 1;

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
profile.bart_pics_cmd       = 'pics  -R T:1024:0:0.1 -R T:2048:0:0.1 -i 20 -S -d5'; %PROUD CINE test-retest settings
    % pics [-l ...] [-r f] [-c] [-s f] [-i d] [-t <string>] [-n] [-g] [-p <string>] [-I ...] [-b d] [-e] [-W <string>] [-d d] [-u f] [-C d] [-f f] [-m ...] [-w f] [-S] [-B d] [-K] <kspace> <sensitivities> <output>

profile.bart_sensemap_cmd   = 'ecalib -I -r20';   
    % ecalib [-t f] [-c f] [-k ...] [-r ...] [-m d] [-S] [-W] [-I] [-1] [-v f] [-a]
    % caldir cal_size  

profile.bart_cc_cmd         = 'cc -p 8 -E';
    % cc [-p d] [-M] [-r ...] [-A] [-S ...] [-G ...] [-E ...] <kspace> <coeff>|<proj_kspace>
    % Types - S: SVD; G: Geometric ;E: ESPIRiT
    
profile.bart_parpool             = 4;
    % 0:     NO  parallel computing in CS reconstruction
    % N > 0: Use parallel computing in CS reconstruction with N workers

%% % CARDIAC PARAMETER
profile.card_RetroHoleInterpolationCenter = 'No'; % 'Yes', 'No'
profile.Cardiac_RetroHoleInterpolation  = 'No'; % 'No', 'Nearest', 'Average', 'Linear', 'Cubic'  (If yes, select 'Average' for best result, is what Lukas said) 
profile.Cardiac_Synchronization         = 'Retrospective'; % 'Retrospective' , 'None'
profile.Cardiac_RetroBinning            = 'Relative'; % 'Relative', 'Absolute', 'None'
profile.Cardiac_RetroPhases             = 15; % for relative binning
profile.Cardiac_HeartPhaseInterval      = 40; % for absolute binning ?
% profile.Cardiac_RetroEndSystoleMs       = 350;
% profile.Cardiac_PhaseWindow             = '';
% profile.Cardiac_RespSync                = 'No';
% profile.Cardiac_RespComp                = 'No';
% profile.Cardiac_RNAV                    = '';

%% % RESPIRATORY PARAMETER
profile.resp_RetroRespiratoryBinning        = 'Yes';  %'Yes', 'No'  
profile.resp_RetroRespiratoryBinningMethod  = 'SG' ;   % 'SG', 'VitalEye'
profile.resp_Phases                         = 4;  %number of respiratory bins

%% % RECON PARAMETER
profile.Recon_CoilCombination               = 'pc';
profile.Recon_ImageSpaceZeroFill            = 'Yes';
profile.Recon_kSpaceZeroFill                = 'Yes';
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
%% % ENCODING PARAMETERS
% profile.Encoding_KzOversampling             = 1; % 

%% % Retrospective Undersampling
% profile.RU_do   = 0; % 1: do retrospective undersampling
% profile.RU_mask = ''; % file location needed, if RU_do = 1

%% save profile
save_profile(profile);
