function profile = load_profile( filename )
%LOAD_PROFILE Summary of this function goes here
%   Detailed explanation goes here

[PATH,~,ext] = fileparts(filename);

if strcmp(PATH, '')
    filename = fullfile(fileparts(fileparts(mfilename('fullpath'))),'profileQueue',filename); 
end

% make script unaffected by USER
C=strsplit(filename,'/');
C{3} = getenv('USER');
filename = strjoin(C,'/');

if strcmp(ext, '.mat')
	profile = load(filename);
elseif any(strcmp(ext, {'.txt' '.dat' '.csv'})),
    % re-format the cell-table into struct
	%T = readtable(filename, 'Delimiter','\t','Format','auto'); renske
    T = readtable(filename, 'Delimiter','\t','Format','auto');
	profile = {};
    varnames = T.Properties.VariableNames;   
    for ii=1:size(T,1)
        if eval(sprintf('isnan(str2double(T.%s{ii}))',varnames{2}))
            eval(sprintf('profile.(char(T.field{ii})) = T.%s{ii};',varnames{2})) 
        elseif eval(sprintf('ischar(T.%s{ii})',varnames{2}))
            eval(sprintf('profile.(char(T.field{ii})) = str2double(T.%s{ii}); ',varnames{2}))
        else
            eval(sprintf('profile.(char(T.field{ii})) = T.%s{ii}; ',varnames{2}))
        end
    end
    for jj=1:size(T,2)-1
        if eval(sprintf('isa(T.%s,''double'')',varnames{jj+1}))
            eval(sprintf('profile.(char(T.field{ii})) = [profile.(char(T.field{ii})) T.%s(ii)]; ',varnames{jj+1}))
        end
    end
end

end