classdef Logger < handle
	%LOG Summary of this class goes here
	%   Detailed explanation goes here
	
	properties (Constant = true)
		fmt = '%08.3f s';
        fmt2 = '%02ih : %02im : %02is';
		width = 48;
	end
	
	properties (SetAccess = immutable)
		quiet = false;
		t0 = 0;
	end
	
	properties (Hidden = true)
		t = 0;
		line = '';
	end
	
	properties
		text = {};
        savename = datestr(datetime('now'),30);
	end
	
	methods
		function self = Logger(mode)
			if nargin < 1; elseif strcmp(mode, 'quiet')
				self.quiet = true;
            end
            % diary test
            PATH = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'logs');
            if ~isdir(PATH), mkdir(PATH), end
            diary(fullfile(PATH,[self.savename,'_diary.log']));
			disp('[00h : 00m : 00s]  START CS reconstruction:');
            self.text{1} = '[00h : 00m : 00s]  START CS reconstruction:';
			self.t0 = tic;
		end
		
		function start(self, str, savename)
			self.t = toc(self.t0);
            self.savename = savename;
% 			time = sprintf(['[' self.fmt ']'], self.t);
            time = sprintf(['[' self.fmt2 ']'], gettime(self.t));
			self.line = [time '  ' str];
			self.text{end+1} = self.line;
			if ~self.quiet; disp(self.line); end
            save(self)
		end
		
		function stop(self)
% 			time = sprintf(['[' self.fmt ']'], toc(self.t0)-self.t);
            time = sprintf(['[' self.fmt2 ']'], gettime(toc(self.t0)-self.t));
            padding = repmat(' ', [1 self.width]);
			self.line = [padding ' ' time];
			self.text{end+1} = self.line;
			if ~self.quiet; disp(self.line); end
            save(self)
		end
		
		function note(self, val, txt)
%  			time = sprintf(['[' self.fmt ']'], toc(self.t0));
			time = sprintf(['[' self.fmt2 ']'], gettime(toc(self.t0)));
            if nargin < 2+1; txt = ''; end
			str = deblank(evalc('disp(val)'));
			self.line = sprintf('%s %s', [time ' ' txt], str);
			self.text{end+1} = self.line;
			if ~self.quiet; disp(self.line); end
            save(self)
		end
		
		function finish(self)
% 			time = sprintf(['[' self.fmt ']'], toc(self.t0));
            time = sprintf(['[' self.fmt2 ']'], gettime(toc(self.t0)));
            self.text{end+1} = [time ' END CS reconstruction.'];
			if ~self.quiet; disp(time); end
            save(self)
			disp('END reconstruction.');
            % diary test
            diary('off');
		end
		
		function print(self)
			for ii = 1:length(self.text)
				disp(self.text{ii});
			end
		end
		
		function save(self)
			PATH = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'logs');
			file = fopen(fullfile(PATH, [self.savename, '.log']),'w');
			for ii = 1:length(self.text);
				fprintf(file,'%s\n',self.text{ii});
			end
			fclose(file);
		end
		
	end
	
end

function out = gettime(t)
    hh = floor(t / 3600);
    mm = floor(t/60)-hh*60;
    ss = floor(t-hh*3600-mm*60);
    out = [hh, mm, ss];
end