function DIR = dir_out( RP )
%DIR_OUT Summary of this function goes here
%   Detailed explanation goes here

% v.q.pronk@amc.uva.nl; initial code
% edited by: l.m.gottwald@amc.nl


DIR = fullfile(RP.data_dir, ['recon_out_',RP.name]);

if ~exist(DIR, 'dir')
  mkdir(DIR);
end

end

