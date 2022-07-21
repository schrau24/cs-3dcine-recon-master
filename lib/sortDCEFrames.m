function [extr2, nt] = sortDCEFrames(mrecon, DCETimeResolution)
% build profile timing table for sequence with SPAIR/SPIR pulses

noisemeas = sum(mrecon.Parameter.Labels.Index.typ == 5);
phasecorrmeas = sum(mrecon.Parameter.Labels.Index.typ == 3);
data_in_use_index = false(size(mrecon.Parameter.Labels.Index.rtop));
data_in_use_index(1+noisemeas+phasecorrmeas:mrecon.Parameter.Labels.OriginalLabelLength)=true;  %rm: you already only read in data of typ ==1?

dataLength = size(mrecon.Data,2);
% channels
chan=mrecon.Parameter.Labels.CoilNrs(:,1);
nchan = length(chan);

scanDuration = mrecon.Parameter.Labels.ScanDuration;
% if fat sat used, we need to calculate true readout times throughout the
% scan. otherwise, simply use TR and dataLength
disp(['Sorting data into ' num2str(DCETimeResolution) ' s frames'])
if ~mrecon.Parameter.Labels.SPIR
    % the TR and expected number of readouts per DCE frame
    TR = mrecon.Parameter.Labels.RepetitionTime/1000;                   % in seconds
    nROperFrame = floor(DCETimeResolution/TR);
else
    
    % need to get very specific with timing:
    % what is the dummy duration (post base and pre SPAIR timing)
    tmp = mrecon.Parameter.GetObject('SQ`dummy');
    dummyDur = tmp.GetValue('dur')/1000;
    % and what is the SPAIR duration
    tmp = mrecon.Parameter.GetObject('SQ`spir');
    SPAIRDur = tmp.GetValue('dur')/1000;
    
    TR = mrecon.Parameter.Labels.RepetitionTime/1000;
    
    % TFE parameters, nost unused
    tfe_shot_dur = mrecon.Parameter.GetValue('IF_tfe_shot_dur')/1000;   % in seconds
    SPAIR_tr = mrecon.Parameter.GetValue('IF_SPIR_adiab_spir_rep_time')/1000;
    tfe_shots = mrecon.Parameter.GetValue('IF_tfe_shots');
    
    % data is only acquired during the last portion of the tfe_shot
    tfe_acq_dur = mrecon.Parameter.GetValue('IF_tfe_acq_dur')/1000;     % in seconds
    tfe_factor = mrecon.Parameter.GetValue('IF_tfe_factor');
    
    % now we can build up the timing table
    currTime = 0;   % tracking time, in s
    timingTable = zeros(tfe_shots*tfe_factor,2);
    for shot = 1:tfe_shots
        if shot == 1
            % note there are 6! dummy+SPAIR pulses before data acquisition begins
            % more accurately, with calibration scans, the first profile
            % acquisition is not until 1951.52 ms
            %             currTime = currTime + 6*(dummyDur+SPAIRDur);
            currTime = currTime + 1.95152;
        end
        
        % the indices for profiles of the current shot
        idx = (1:tfe_factor)+((shot-1)*tfe_factor);
        
        timingTable(idx,1) = idx;
        timingTable(idx,2) = currTime+(1:tfe_factor)'*TR;
        
        % update currTime
        currTime = currTime + tfe_factor*TR + (dummyDur+SPAIRDur);
    end
    
    % assign frames based on timings
    diffTimingTable = timingTable(:,2) - timingTable(1,2);
    nt = floor(max(diffTimingTable/DCETimeResolution));
    
    clear sortedFrames
    for frame = 1:nt
        if frame == 1
            idx = find(diffTimingTable < DCETimeResolution);
        elseif frame < nt
            idx = find(diffTimingTable >= DCETimeResolution*(frame-1) & ...
                diffTimingTable < DCETimeResolution*(frame));
        else
            idx = find(diffTimingTable >= DCETimeResolution*(frame-1));
        end
        sortedFrames{frame} = idx;
    end
end
disp(['Resulting in ' num2str(nt) ' frames'])

% WRITE LABEL
extr2 = zeros(size(data_in_use_index));
startIdx = noisemeas + phasecorrmeas + 1;
for ii=1:length(sortedFrames)
    if ii==1
        idx = 1:(max(sortedFrames{ii})*nchan);
    elseif ii < nt
        idx = (max(sortedFrames{ii-1})*nchan+1):(max(sortedFrames{ii})*nchan);
    else
        idx = (max(sortedFrames{ii-1})*nchan+1):length(extr2);
    end
    idx = idx+startIdx;
    extr2(idx) = ii-1;
end