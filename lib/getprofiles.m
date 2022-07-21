function filenames = getprofiles()

EXTENSIONS = {
	'*.pro.dat', 'Reconstruction profile human-readable (.pro.dat)';
	'*.pro.mat', 'Reconstruction profile binary  (.pro.mat)';
};

TITLE = 'Select Reconstruction Profile(s)';
PATH = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'profileQueue');

filenames = uigetfile(EXTENSIONS, TITLE, PATH, 'MultiSelect', 'on');
if ischar(filenames), filenames = {filenames}; end;

end
