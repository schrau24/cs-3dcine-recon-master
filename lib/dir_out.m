function DIR = dir_out( RP )
%DIR_OUT Summary of this function goes here
%   Detailed explanation goes here

% v.q.pronk@amc.uva.nl; initial code
% edited by: l.m.gottwald@amc.nl
if isfield(RP, 'output_dir')
    outputdir = RP.output_dir;
else
    outputdir = RP.data_dir;
end

DIR = fullfile(outputdir, ['recon_out_',RP.name]);

if ~exist(DIR, 'dir')
  mkdir(DIR);
end

end

