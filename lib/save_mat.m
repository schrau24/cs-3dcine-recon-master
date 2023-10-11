function save_mat( path, filename, variable, value, mode )
%SAVE_MAT Summary of this function goes here
%   Detailed explanation goes here

if ~exist('mode', 'var'), mode = 'v7.3'; end;

str = struct();
str.(variable) = value;
save(fullfile(path, [filename '.mat']), '-struct', 'str', ['-', mode]);

end
