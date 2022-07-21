function [rnav] = RNAVExtraction(labels)

%% find RNAV in MRecon.Data, copied from add_rnav_to_labels.m
if isempty(labels.RNAV)
    error('ERROR in add_rnav_to_labels: rnav array is empty');
end

noisemeas = sum(labels.Index.typ == 5);
phasecorrmeas = sum(labels.Index.typ == 3);

data_in_use_index = false(size(labels.Index.rtop));
data_in_use_index(1+noisemeas+phasecorrmeas:labels.OriginalLabelLength)=true;

rtops = labels.Index.rtop(data_in_use_index); 
rr = labels.Index.rr(data_in_use_index); 
rtops_diff = zeros(size(rtops));
rtops_diff(2:end) = double(rtops(2:end))-double(rtops(1:end-1));
neg_indx = rtops_diff<0;
rtops_diff(neg_indx) = rtops_diff(neg_indx)+double(rr(neg_indx([2:end,1],1)));

rnavs = (rtops_diff>20)&(rtops_diff<50);% RNAV takes ~40ms
rnavs_labelnr = find(rnavs);
beams_all = [labels.RNAV.Beam1.PrepLead,labels.RNAV.Beam1.AcqLead]';
beams_all(abs(beams_all)>100)=0;
beams_all=beams_all(find(beams_all,1,'first'):find(beams_all,1,'last')); % remove trailing zeros
beams_all(1:end-size(rnavs_labelnr,1))=[]; % LG: AcqLead must match length of rnavs_labelnr

% create rnav numbers for index list
rnav_number = zeros(size(rtops_diff));
rnav_beamloc = zeros(size(rtops_diff));
for ii=1:size(rnavs_labelnr)
    if ii==1
        rnav_number(1:rnavs_labelnr(ii))=ii;
        rnav_beamloc(1:rnavs_labelnr(ii))=beams_all(ii);
    else
        rnav_number(rnavs_labelnr(ii-1)+1:rnavs_labelnr(ii))=ii;
        rnav_beamloc(rnavs_labelnr(ii-1)+1:rnavs_labelnr(ii))=beams_all(ii);
    end
end

% smooth the resp nav signal for better delineation of resp phases
% later, use half the number of points sampled between each rnav
rnav_beamloc = smooth(rnav_beamloc,round(mean(diff(rnavs_labelnr))/2));

% get rid of errant high rnav signals
rnav_beamloc(rnav_beamloc > 25) = 0;

% % normalize to get end-expiration at 0 by subtracting the mean of the peaks
% [~,idx] = findpeaks(rnav_beamloc);
% rnav_beamloc = rnav_beamloc - mean(rnav_beamloc(idx));

rnav = zeros(size(labels.Index.rtop));
rnav(data_in_use_index) = rnav_beamloc;