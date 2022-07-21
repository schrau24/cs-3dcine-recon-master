function move_pro( filename, profile, mode )
%MOVEPROFILE Summary of this function goes here
%   Detailed explanation goes here

from = fullfile(fileparts(fileparts(mfilename('fullpath'))),'profileQueue', filename);
to = fullfile(dir_out(profile), filename);
archivedir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'profileQueue','Archive');
if ~exist(archivedir, 'dir')
  mkdir(archivedir);
end
archivedir = fullfile( archivedir, filename);

if exist('mode','var') && strcmp(mode, 'copy')
	copyfile(from, to);
elseif exist('mode','var') && strcmp(mode, 'archive')
    copyfile(from, to);
    movefile(from, archivedir);
else
	movefile(from, to);
end

end
