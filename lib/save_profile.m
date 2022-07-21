function save_profile( profile, ext )
%SAVE_PRO Summary of this function goes here
%   Detailed explanation goes here

if ~exist('ext','var'), ext = '.dat'; end;
folderpath = fullfile(fileparts(fileparts(mfilename('fullpath'))),'profileQueue');
if ~isdir(folderpath) , mkdir(folderpath); end;
file = fullfile(folderpath, [profile.name '.pro' ext]);

if strcmp(ext, '.mat')
	save(file,'-struct','profile');
elseif any(strcmp(ext, {'.txt' '.dat' '.csv'}))
	S = profile;
	K = fieldnames(S);
	V = cell(length(S)); for ii=1:length(K), V{ii}=S.(K{ii}); end;
	T = table(K, V.', 'VariableNames', {'field' 'value'});
	
	writetable(T, file, 'Delimiter', '\t');
else
	error('The extension "%s" is not a valid format', ext);
end

end