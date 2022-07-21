function [extr2, nt] = sortDCEFrames(mrecon, DCETimeResolution)

noisemeas = sum(mrecon.Parameter.Labels.Index.typ == 5);
phasecorrmeas = sum(mrecon.Parameter.Labels.Index.typ == 3);
data_in_use_index = false(size(mrecon.Parameter.Labels.Index.rtop));
data_in_use_index(1+noisemeas+phasecorrmeas:mrecon.Parameter.Labels.OriginalLabelLength)=true;  %rm: you already only read in data of typ ==1?

dataLength = size(mrecon.Data,2);
% channels
chan=mrecon.Parameter.Labels.CoilNrs(:,1);
nchan = length(chan);

% if fat sat used, we need to calculate true readout times throughout the
% scan. otherwise, simply use TR and dataLength
disp(['Sorting data into ' num2str(DCETimeResolution) ' s frames'])
if ~mrecon.Parameter.Labels.SPIR
    % the TR and expected number of readouts per DCE frame
    TR = mrecon.Parameter.Labels.RepetitionTime/1000;                   % in seconds
    nROperFrame = floor(DCETimeResolution/TR);
else
    % TFE parameters
    tfe_shot_dur = mrecon.Parameter.GetValue('IF_tfe_shot_dur')/1000;   % in seconds
    tfe_shots = mrecon.Parameter.GetValue('IF_tfe_shots');
    
    % data is only acquired during the last portion of the tfe_shot
    tfe_acq_dur = mrecon.Parameter.GetValue('IF_tfe_acq_dur')/1000;     % in seconds
    tfe_factor = mrecon.Parameter.GetValue('IF_tfe_factor');
    
    nShotsperFrame = ceil(DCETimeResolution/tfe_shot_dur);
    trueDCETimeResolution = nShotsperFrame*tfe_shot_dur;
    nROperFrame = nShotsperFrame*tfe_factor;
    disp(['With TFE factor = ' num2str(tfe_factor)])
    disp(['True frame length = ' num2str(trueDCETimeResolution) ' s'])
end
% calculate the number of frames in our data
nt = floor(dataLength/nchan/nROperFrame);
disp(['Resulting in ' num2str(nt) ' frames'])

% WRITE LABEL
extr2 = zeros(size(data_in_use_index));
startIdx = noisemeas + phasecorrmeas + 1;
for ii=1:nt
    if ii==1
        idx = startIdx:(nROperFrame*nchan);
    elseif ii < nt
        idx = (nROperFrame*nchan*(ii-1)+startIdx):(nROperFrame*nchan*ii);
    else
        idx = (nROperFrame*nchan*(ii-1)+startIdx):length(extr2);
    end
    extr2(idx) = ii-1;
end