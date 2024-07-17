%% CS-4Dflow-recon
% Run script (F5) to start reconstructions

% Editors (alphabetical order):
% l.m.gottwald@amc.uva.nl
% e.s.peper@amc.uva.nl
% v.q.pronk

%% housekeeping / initialization
restoredefaultpath;
addpath(genpath(fullfile(cd,'./lib/')));
clear;
close all;
clc;

%% run recons
global logger

PROFILES = dir('profileQueue/*.pro.dat');
PROFILES = {PROFILES.name};

for ii = 1:length(PROFILES)
    TIME = datestr(datetime('now'),30);
    logger = Logger();
    try
        pfile = PROFILES{ii};
        profile = load_profile(pfile);
        % update profile name to include today's date
        d = num2str(yyyymmdd(datetime));
        
        if strmatch(profile.resp_ReconType,'resp_motionCorrection')
            profile.name = [profile.name '_resp_motionCorrection' '_' d];
            disp('Running reconstruction:')
            disp([profile.name])
            
            logger.start(profile.name,TIME);
            % Respiratory recon and registration to get offsets
            recon_RespRegister( profile )
            recon_WH_3DCINE_wRespCorrection( profile );
        else
            
            profile.name = [profile.name '_resp_5DRecon' '_' d];
            disp('Running reconstruction:')
            disp([profile.name])
            
            logger.start(profile.name,TIME);
            Recon_5D(profile);
        end 
        
%         profile.name = [profile.name 'resp_5DRecon_plusMotionCorrection' '_' d];
%         disp('Running reconstruction:')
%         disp([profile.name])
%         
%         logger.start(profile.name,TIME);
%         % Respiratory recon and registration to get offsets
%         recon_RespRegister( profile )
%         recon_5D_wTranslations(profile);
        move_pro(pfile, profile, 'archive');
    catch ME
        logger.note(sprintf('WARNING: Error occured in reconstruction! Recon might be incomplete or incorrect!\n%s',ME.message'));
    end
    logger.finish();
    try
%         sendmail_from_amc_textfile('e.m.schrauben@amsterdamumc.nl@amc.uva.nl',...
%             sprintf('Finished recon: %s',pfile),...
%             fullfile(cd,'logs',[logger.savename,'.log']));
    catch ME
        logger.note(sprintf('WARNING: Error occured in sendmail. Log is not sent.\n%s',ME.message'));
    end
end