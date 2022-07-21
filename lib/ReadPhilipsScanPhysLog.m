
function [ DATA HDR ] = ReadPhilipsScanPhysLog(filename, channels, skipprep)
%READPHILIPSSCANPHYSLOG load Philips MRI physiolog file
%
%   [ DATA HDR ] = READPHILIPSSCANPHYSLOG(FILENAME)
%   Reads the sample information (all 10 channels) from FILENAME. The 
%   samples are loaded as signed 32-bit integer values (int32). Sampling
%   rate is 496.03Hz for wireless VCG systems and 500Hz for the older wired sytems.
%
%   [ DATA HDR ] = READPHILIPSSCANPHYSLOG(FILENAME, CHANNELS)
%   Also read the log data, but only loads the signals specified by cell
%   structure CHANNELS. {v1raw, v2raw, v1, v2, ppu, resp, gx, gy, gz, mark}
%   Use {} or 'all' to load all channels (default), or an arbitrary invalid
%   name (i.e. 'none') to ignore the channel samples and only load markers.
%
%   [ DATA HDR ] = READPHILIPSSCANPHYSLOG(FILENAME, CHANNELS, SKIPPREP)
%   Also read the log data and specified channels, but skip preperation
%   fase (marked by empty comment (#) line)
%
%   DATA    is a structure with elements for the channel matrix C, 
%           marker table M and event index vectors I. 
%   DATA.C  32-bit integer matrix containing column vectors for
%           all resp. channels values. (i.e. each signal channel is stored 
%           as column vector). Note that the mark column contains
%           hexadecimal values (this will always be the last column).
%   DATA.M  is a 2-column table containing {marker,index} pairs of the
%           the markers detected in the marker channel. The first
%           column of M contains the mark values >0, and the second column 
%           contains the sample indices of the corresponding mark-value 
%           in the first column. Markers are bit-encoded values (i.e. each
%           event type corresponds to a single bit position). Currently the
%           following bits are known:
%                mark      bit pattern         description
%               0x0001 = 0000.0000.0000.0001 = Trigger ECG
%               0x0002 = 0000.0000.0000.0010 = Trigger PPU
%               0x0004 = 0000.0000.0000.0100 = Trigger Respiration
%               0x0008 = 0000.0000.0000.1000 = Measurement ('slice onset')
%               0x0010 = 0000.0000.0001.0000 = start of scan sequence (decimal 16)
%               0x0020 = 0000.0000.0010.0000 = end of scan sequence (decimal 32)
%               0x0040 = 0000.0000.0100.0000 = Trigger external
%               0x0080 = 0000.0000.1000.0000 = Calibration
%               0x0100 = 0000.0001.0000.0000 = Manual start
%               0x8000 = 1000.0000.0000.0000 = Reference ECG Trigger
%
%           Note that events can occur at the same sample, so bits can be
%           combined and give marker values like 3 (=1+2). The events in
%           DATA.I are extracted from DATA.M by bit masking, so all event
%           index vectors include overlapping events properly.
%
%           Note2: The availablility of markers depend on the software release of the scanner.
%
%   DATA.I  a structure holding index vectors for all known event types.
%           The event indices are derived by scanning the bit-pattern
%           in the marker column.
%   DATA.I.VcgOnset     vector with sample indices of gated events
%   DATA.I.PpuOnset     vector with sample indices of ... events
%   DATA.I.Measurement  vector with sample indices of slice onset events
%   DATA.I.TriggerResp  vector with sample indices of respiration trigger
%   DATA.I.ScannerStart vector with sample indices of start of MRI sequence (incl. preparation)
%   DATA.I.ScannerStop  vector with sample indices of end of MRI sequence
%   DATA.I.TriggerExt   vector with sample indices of external triggers
%   DATA.I.Calibration  vector with sample indices of calibration
%   DATA.I.RefTriggerVcg vector with sample indices of ecg reference trigger
%   HDR     Optional return value for a structure containing all
%           information in the file header. 
%   HDR.ID
%   HDR.ID.Site         site name (from first line)
%   HDR.ID.Release      software release name
%   HDR.ID.SWID         software number
%   HDR.DATETIME        year,month,date,hour,min,sec
%   HDR.STATS           statistics from line 3
%   HDR.DockableTable   true/false
%   HDR.COLUMN_NAMES    resp. signal names: {v1raw, v2raw, v1, v2, ppu, 
%                                           resp, gx, gy, gz, mark}
%
%   Examples:
%
%   D = ReadPhilipsScanPhysLog('SCANPHYSLOG20110420142406.log');
%   Reads all physiological channels and marker flags. The marker index
%   information will automatically put into D.M and D.I.
%
%   D = ReadPhilipsScanPhysLog('SCANPHYSLOG20110420142406.log',{'v1','v2','resp'});
%   Reads both artefact corrected VCG and respiration channels.  The marker index
%   information will automatically put into D.M and D.I.
%
%   D = ReadPhilipsScanPhysLog('SCANPHYSLOG20110420142406.log','none');
%   Ignore channels, and only load markers (i.e. DATA.C==[])
%
%   [D H] = ReadPhilipsScanPhysLog('SCANPHYSLOG20110420142406.log');
%   Also parses the header and stores info in H.
%
%   Using start/stop markers to extract epochs:
%
%       beginindex = D.M(find(D.M(:,1)==16),2);
%       endindex   = D.M(find(D.M(:,1)==32),2);
%       epoch = D.C(beginindex:endindex,:);
%
%   However, the start marker is not properly synchronized with the onset of the 
%   first volume. As a workaround you could use the end-marker and real duration of 
%   the scan and calculate the start index yourself:
%
%         freq=496; % for wireless (wBTU) VCG ssytems, 500 for wired
%         TR=2.3;
%         nrvolumes=200;
%         endindex = D.M(find(D.M(:,1)==16),2);
%         beginindex = endindex  - nrvolumes * TR * freq;
%         epoch = logdata.C(beginindex:endindex,:);
%
%   Compatability
%   Developed and tested with Matlab R2009a (7.8.0) - Windows 64-bit.
%   bitand behaviour was changed in Version 7 (R14), so you might have
%   compatibility issues with release before R14.
%
% Copyright 2011 Academical Medical Center Amsterdam
% Created by Paul F.C. Groot
% $Rev::            $:  Revision of last commit 
% $Author::        $:  Author of last commit 
% $Date::             $:  Date of last commit 
%%
% Additonal info from:
% http://dbic.dartmouth.edu/wiki/index.php/Recording_Physiology
% The onset of the file corresponds to 15 seconds prior to the start of the
% scan. There is a single line containing a hash mark (#) and no numbers
% that serves as a time-stamp of the onset of the actual scan. The last
% line of the file corresponds to the exact end of the scan. Note: For
% functional scans, the # mark corresponds to the onset of the template,
% not the onset of the functional run itself, and recording continues
% during any amount of time that passes between the end of the template
% scan and the onset of the actual functional run (which is variable).
% Therefore, it is advisable to time-lock physiological response data
% during functional scans based on when the scan ENDS.
 
    %% Initialize
    if nargin<2
        channels = {};
    elseif ~iscell(channels)
        channels = { channels };
    end
    if nargin<3 || isempty(skipprep)
        skipprep = false;
    end
    
    bParseHeader = nargout>=2;
       
    % open philips log file, and parse all lines
    fid=fopen(filename,'rt');
    if fid==-1
        error('Couldn''t open %s',filename);
    end
    
    %% fill HDR (read header)
    % the file starts with a header, like this one:
    % ## AMC Amsterdam, Release hestia (SWID 118)
    % ## Wed 20-04-2011 14:24:06
    % ## 179 -1351 306 171 -567 266 -30 375 0
    % ## Dockable table = FALSE
    % # v1raw v2raw  v1 v2  ppu resp  gx gy gz mark
	%
	% NB. Release 5 files start with:
	% ## A.M.C AMSTERDAM (SRN = 42151), Release ingeniamn (SWID 31)
	%
    % define regular expressions for the resp. header lines
    if bParseHeader
        p{1} = '##\s*(.*),\sRelease\s(\w+)\s\(SWID (\d+)\)$';
        p{2} = '##\s*[\s\w]*(\d\d)-(\d\d)-(\d\d\d\d)\s(\d\d):(\d\d):(\d\d)$';
        p{3} = '##\s*([\s\-\d])+$';
        p{4} = '##\s*Dockable table\s*=\s*(\w+)$';
    %   p{5} = '#\s*(\w+)...'; See note below '#'
        iPattern = 1;
    end
    
    H = cell(4,1);      % this will hold the results of regexp parsing
    COLUMN_NAMES = [];  % signal name from last header line
    % read all header lines (starting with ## or #) until single #
    while ~feof(fid)
        str = fgetl(fid);
        if strncmp(str,'##',2)
            if bParseHeader
                T = regexp(str, p{iPattern}, 'tokens');
                if length(T)~=1 
                    warning('SCANPHYSLOG:invalidHeader','Invalid header at line 1: %s', str); 
                elseif iPattern<=length(p) 
                    H{iPattern}=T{1}; 
                    iPattern=iPattern+1; % continue with next pattern if current one matched
                else
                    warning('SCANPHYSLOG:invalidHeader','Ignoring additional header line: %s', str);
                end;
            end
        elseif strncmp(str,'#',1)
            % this last header line is parsed a bit differently: this way we can
            % continue without error if the file has a different number of
            % channels (butthe marker should always be the last...)
            COLUMN_NAMES = regexp(str(2:end),'\w+','match'); 
            break;
        else
            error('Invalid header');
        end
    end
    % assemble the returned HDR structure when requested
    if bParseHeader
        % translate file identification info from first header line into struct
        if length(H)>=1 && ~isempty(H{1})
            T = H{1};
            ID = struct('Site', T{1}, 'Release', T{2}, 'SWID', str2double(T{3}) );
        else
            ID = struct('Site', 'Unknown site', 'Release', 'unknown release', 'SWID', 0 );
        end
        % translate scan data/time info from second header line into struct
        if length(H)>=2 && ~isempty(H{2})
            T = H{2};
            DATETIME = struct(  'year',     uint16(str2double(T{3})), ... 
                                'month',    uint16(str2double(T{2})), ...
                                'day',      uint16(str2double(T{1})), ...
                                'hour',     uint16(str2double(T{4})), ... 
                                'min',      uint16(str2double(T{5})), ...
                                'sec',      uint16(str2double(T{6})) ...                                 
                            );
        else
            DATETIME = struct(  'year',     uint16(1900), ... 
                                'month',    uint16(1), ...
                                'day',      uint16(1), ...
                                'hour',     uint16(0), ... 
                                'min',      uint16(0), ...
                                'sec',      uint16(0) ...                                 
                            );
        end
        % translate signals statistics from thirth header line into struct
        if length(H)>=3 && ~isempty(H{3})
            T = H{3};
            STATS = sscanf(T{1},'%d')';
        else
            STATS = NaN .* zeros(1,9);
        end
        % translate DockableTable flag from fourth header line into struct
        if length(H)>=4 && ~isempty(H{4})
            T = H{4};
            DockableTable = strcmpi(T,'TRUE');
        else
            DockableTable = false;
        end
        HDR = struct('ID', ID, 'DATETIME', DATETIME, 'STATS', STATS, ...
            'DockableTable', DockableTable, 'COLUMN_NAMES', COLUMN_NAMES);
    end
    %% skip preparation?
    if skipprep
        startpos = ftell(fid); % remember start position in ase we didn't find a #
        while ~feof(fid)
            str = fgetl(fid);
            if strcmp(strtrim(str),'#')
                break;
            end
        end
        if feof(fid)
            warning('SCANPHYSLOG:prepNotFound','Preparation fase not found: reading all samples');
            fseek(fid,startpos,'bof'); % go back to where we came from
        end
    end
 
    %% fill DATA by loading requested columns
    % Create a logical vector S which is true for selected channels. This
    % vector is needed below with textscan to skip unrequired columns.
    % First make sure that the channel name cell array is properly filled.
    nColumns = length(COLUMN_NAMES);
    if isempty(channels) || sum(strcmpi(channels,'all'))>0
        channels = COLUMN_NAMES; % just include all channels
    end
    % Also create an ordering vector to reshuffle the columns back into the
    % order the caller requested (used below to reshuffle C)
    channel_ordering = []; 
    S = logical(zeros(1,nColumns)); %#ok<LOGL>
    for iChannel=1:length(channels)
        which_channel_vector = strcmpi( channels{iChannel}, COLUMN_NAMES );
        S = S | which_channel_vector;
        iChannelPos = find(which_channel_vector);
        if (iChannelPos>0)
            channel_ordering(end+1) = iChannelPos; %#ok<AGROW>
        end
    end
    % The channel_ordering vector should be renumbered because columns
    % might have been skipped (i.e. [2 5 3] should become [1 3 2].
    % So, loop and replace the lowest by incrementing index.
    column_ordering = zeros(1,length(channel_ordering));
    temp_channel_ordering = channel_ordering;
    for iChannel=1:length(temp_channel_ordering)
        iMin = find(temp_channel_ordering==min(temp_channel_ordering));
        column_ordering(iMin) = iChannel;
        % set to huge value so it can't be the lowest anymore in next
        % iteration:
        temp_channel_ordering(iMin) = Inf; 
    end
 
    % prepare a format string for parsing lines:
    % - use %* to skip columns
    % - use sprintf generate the format string: %c is placeholder for * or null, %% is just %
    % - the last column is a hexadecimal value (which is not supported yet by textscan with %x),
    %   so we first read it as plain string and then convert explicitly.
    %   The bit pattern must be preserved because different event types can occur at the same sample.
    format = [ sprintf('%%%cd ',(S(1:end-1)==0).*'*') '%4s']; % e.g. '%*d %*d %d %*d %d %d %d %d %d %4s
    format = regexprep(format,'%[^*]d','%d'); % remove null characters
    
    % read all sample rows, skipping stray '#' lines along the way
    C=textscan(fid,format,'MultipleDelimsAsOne',1,'CommentStyle','#');
    fclose(fid);
    
    %% translate marker column to event indices (sample indices)
    % get last (mark) column and translate 4 hex digits to unsigned number
    markers = uint16(hex2dec(C{end}));
%   fprintf('markers: '); fprintf(' %d',unique(markers)); fprintf('\n');
    % then mask all known event types (marked by single bit) and create
    % index vectors for each type; create a struct along the way
    n = length(markers);
    onevec = ones(n,1,'uint16');
    I.VcgOnset      = uint32(find(bitand(markers,1*onevec))); % VCG wave onset
    I.PpuOnset      = uint32(find(bitand(markers,2*onevec))); % PPU wave onset
    I.TriggerResp   = uint32(find(bitand(markers,4*onevec))); % Respiration trigger
    I.Measurement   = uint32(find(bitand(markers,8*onevec)));
    I.ScannerStart  = uint32(find(bitand(markers,16*onevec)));
    I.ScannerStop   = uint32(find(bitand(markers,32*onevec)));
    I.TriggerExt    = uint32(find(bitand(markers,64*onevec)));
    I.Calibration   = uint32(find(bitand(markers,128*onevec)));
    I.RefTriggerVcg = uint32(find(bitand(markers,32768*onevec)));
    % convert marker channel to numerical values if this column should be
    % included in the output 
    if S(end)
        C{end} = int32(markers); % all channels should be the same int32
    else
        C = C(1:end-1); % just remove markers from output
    end
    % create matrix from columns in cell's and reshuffle columns in order caller requested
    C = cell2mat(C(column_ordering));
    % instead of simply returning the marker vector, a table is created
    % which has marker values in the first column and indices (sample
    % numbers) in the second column
    iMarkerIndices = uint32(find(markers>0));
    M = [uint32(markers(iMarkerIndices)), iMarkerIndices];
    
    %% and finally assemble DATA struct to return
    DATA = struct('C', C, 'M', M, 'I', I);
end
