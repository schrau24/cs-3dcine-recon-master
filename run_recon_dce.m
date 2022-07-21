%% CS-recon
% Run script (F5) to start reconstructions

% Editors (alphabetical order):
% l.m.gottwald@amc.uva.nl
% e.s.peper@amc.uva.nl
% v.q.pronk

% e.m.schrauben@amsterdamumc.nl, updated to run respiratory registration
% and correction of coronaries

%% housekeeping / initialization
restoredefaultpath;
addpath(fullfile(cd,'./lib/'));
% addpath(genpath(fullfile(cd,'./Elastix/')));
clear;
close all;
clc;

%% run recons
global logger
PROFILES = getprofiles();
for ii = 1:length(PROFILES)
    TIME = datestr(datetime('now'),30);
    logger = Logger();
    try
        pfile = PROFILES{ii};
        profile = load_profile(pfile);
        
        % update profile name to include today's date
        d = num2str(yyyymmdd(datetime));
        profile.name = [profile.name '_' d];
        disp(['Running reconstruction:'])
        disp([profile.name])
        
        logger.start(profile.name,TIME);

        recon_dce( profile );
        move_pro(pfile, profile, 'archive');
    catch ME
        logger.note(sprintf('WARNING: Error occured in reconstruction! Recon might be incomplete or incorrect!\n%s',ME.message'));
    end
    logger.finish();
    try
%         sendmail_from_amc_textfile('e.m.schrauben@amsterdamumc.nl',...
%             sprintf('Finished recon: %s',pfile),...
%             fullfile(cd,'logs',[logger.savename,'.log']));
    catch ME
        logger.note(sprintf('WARNING: Error occured in sendmail. Log is not sent.\n%s',ME.message'));
    end
end