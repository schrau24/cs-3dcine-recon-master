function pear_signal = readPEARSignal(RESP,RP)
% Grab, load, and cut PEAR signal to match the length of our scan
% More details can be found in ReadPhilipsScanPhysLog.m

% emschrauben@amsterdamumc.nl

% RESP 		- the MRecon object containing acquisition info 
% RP 		- the PROUD profile used to help find the log
% pear_signal 	- the pear signal at each readout in the acquisition

% Have the user select the physlog file, then read in the PEAR signal
[fName, pName] = uigetfile([RP.data_dir '/*.log'],...
    'select the phys log containing PEAR belt');
D = ReadPhilipsScanPhysLog(fullfile(pName,fName),{'resp'},1);

% this is the set frequency for Philips wireles phys monitoring systems
freq = 496.03;
% plot(1/freq*(1:length(D.C)),D.C)

% start time of the PEAR signal does not perfectly align with the true
% acquisition timing, so we need to define it manually here based on the
% TR and number of readouts
chan = RESP.Parameter.Labels.CoilNrs(:,1);
nchan = length(chan);
TR = RESP.Parameter.Scan.TR / 1000;
nReadouts = numel(find(RESP.Parameter.Labels.Index.chan==chan(1)));
endindex = length(D.C);
beginindex = (endindex - floor(nReadouts * TR * freq));
PEARsignal = D.C(beginindex:endindex);

% now we use nearest-neighbor to get resp locations for all of our readouts

% these two timings should be nearly equal
acqTime = 0:TR:(nReadouts-1)*TR;
PEARTime = 1/freq * (0:length(PEARsignal)-1);

resp_PEAR = PEARsignal(knnsearch(PEARTime',acqTime'));

%% WRITE LABEL
noisemeas = sum(RESP.Parameter.Labels.Index.typ == 5);
phasecorrmeas = sum(RESP.Parameter.Labels.Index.typ == 3);
data_in_use_index = false(size(RESP.Parameter.Labels.Index.rtop));
data_in_use_index(1+noisemeas+phasecorrmeas:RESP.Parameter.Labels.OriginalLabelLength)=true;

pear_signal = zeros(size(data_in_use_index)) + mean(resp_PEAR);
for ii=1:length(resp_PEAR)
    %     if ii==1
    %         pear_signal(1:nchan)=resp_PEAR(ii);
    %     elseif ii < length(resp_PEAR)
    ind = ((ii-1)*nchan + 1):(ii*nchan);
    pear_signal(ind)=resp_PEAR(ii);
    %     else
    %         pear_signal(sg_sampling_index(ii):end)=sg_signal_bp_mean(ii,1);
    %     end
end
