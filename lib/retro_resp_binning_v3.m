function labels = retro_resp_binning_v3(labels,method)
% labels = retro_resp_binning_new(labels)
%
% retrospective respiratory binning
% labels: mrecon.Parameter.Labels
% method: method == 0   RNAV binning
% method: method == 1   VitalEye binning
% method: method == 2   self gating
%
% 20-Feb-2020 l.m.gottwald@amsterdamumc.nl

%% get respiratory data
nchan = size(labels.CoilNrsPerStack{1,1},1);
nenc = max(labels.Index.extr1)+1;
noisemeas = sum(labels.Index.typ == 5);
phasecorrmeas = sum(labels.Index.typ == 3);
data_in_use_index(1+noisemeas+phasecorrmeas:labels.OriginalLabelLength)=true;

if method == 0 % RNAV
    resp_signal = double(labels.Index.rnav(data_in_use_index));
    resp_signal = resp_signal - median(resp_signal);
elseif method == 1 % VitalEye
    resp_signal = -double(labels.Index.na(data_in_use_index));
    % fill leading zeros
    startindex = find(resp_signal,1,'first');
    resp_signal(1:startindex)=resp_signal(startindex);    
    resp_signal = resp_signal + 3.3*1e4;
elseif method == 2 % SG
    resp_signal = double(labels.Index.sg(data_in_use_index));    
else
    error('ERROR in retro_resp_binning: no method defined.');
end

%% prepare respiratory signal

% remove unnecessary samples
step = double(nchan*nenc);
resp_signal_steps = [1:step:size(resp_signal,1)]';
resp_signal_short = resp_signal(resp_signal_steps);

% smooth signal over 1 second
time_tr = labels.RepetitionTime;
sma_val = double(round(1000/nenc/time_tr));
resp_signal_short = smooth(resp_signal_short,sma_val);

%% find extrema in signal

% define minimal time distance
max_bpm = 45; % maximal breath per minute
min_seperation = round(60000/max_bpm/nenc/time_tr);

% calculate minimal peak prominence 
% as 1/4 of the median maxima/minima distance
maxima_loc = find(islocalmax(resp_signal_short,'MinSeparation',min_seperation));
minima_loc = find(islocalmax(-resp_signal_short,'MinSeparation',min_seperation));
minmax_length = min(length(maxima_loc),length(minima_loc))-1;
minmaxdiff = resp_signal_short(minima_loc(2:minmax_length))+resp_signal_short(maxima_loc(2:minmax_length));
mpp = median(abs(minmaxdiff))*0.25;

% calculate extrema
maxima_loc = find(islocalmax(resp_signal_short,'MinSeparation',min_seperation,'MinProminence',mpp));
minima_loc = find(islocalmax(-resp_signal_short,'MinSeparation',min_seperation,'MinProminence',mpp));
[minima_loc, maxima_loc] = correct_double_extrema(minima_loc,maxima_loc,resp_signal_short);

% phase bin the respiratory signal in 100 bins 
% (0: peak inspiration, 100: peak expiration)
nbins = 100;
phase_bin = phase_binning(resp_signal_short,maxima_loc, minima_loc, nbins);

resp_signal_typ = false(size(resp_signal_short)); % by default rejected
resp_signal_bin = zeros(size(resp_signal_short)); % by default bin #0

limits=[];
limits(1)=100; % expiration state
limits(2)=40; % 60% acceptance rate / rejected
% limits(3)=30; % 70% acceptance rate / rejected

for ilimit = 1:size(limits,2)-1
    idx = phase_bin<=limits(ilimit) & phase_bin>limits(ilimit+1);
    resp_signal_typ(idx)=true;
    resp_signal_bin(idx)= ilimit;
end


%% %%%%%%
% figure,hold on
% t = 1:length(resp_signal_short);
% plot(resp_signal_short)
% plot(t(maxima_loc),resp_signal_short(maxima_loc),'bo');
% plot(t(minima_loc),resp_signal_short(minima_loc),'ro');

%% Rescale
% allocate temp. arrays
resp_signal_temp=resp_signal;
phase_bin_temp=resp_signal;
resp_signal_bin_temp=resp_signal;
resp_signal_typ_temp=resp_signal;
for ii=1:size(resp_signal_steps,1)
    resp_signal_temp(resp_signal_steps(ii):resp_signal_steps(ii)+step-1)=resp_signal_short((ii));
    phase_bin_temp(resp_signal_steps(ii):resp_signal_steps(ii)+step-1)=phase_bin((ii));
    resp_signal_bin_temp(resp_signal_steps(ii):resp_signal_steps(ii)+step-1)=resp_signal_bin((ii));
    resp_signal_typ_temp(resp_signal_steps(ii):resp_signal_steps(ii)+step-1)=resp_signal_typ((ii));
end
resp_signal = resp_signal_temp;
phase_bin = phase_bin_temp;
resp_signal_bin = resp_signal_bin_temp;
resp_signal_typ = resp_signal_typ_temp;
clear resp_signal_norm_temp phase_bin_temp resp_signal_bin_temp resp_signal_typ_temp

%% Overwrite labels
typ = labels.Index.typ(data_in_use_index);
typ(~resp_signal_typ) = 2; % label data out of respiratory limits as rejected
labels.Index.typ(data_in_use_index) = typ;
labels.Index.extr2(data_in_use_index) = uint16(resp_signal_bin); % LG: use extr2 dimension for resp binning

labels.Index.resp_phase = zeros(size(labels.Index.typ));
labels.Index.resp_phase(data_in_use_index) = phase_bin;
labels.Index.resp_signal = zeros(size(labels.Index.typ));
labels.Index.resp_signal(data_in_use_index) = resp_signal;

%% plot signal
%{
figure
subplot(311), hold on
plot(resp_signal)
subplot(312), hold on
t = 1:length(resp_signal_short);
plot(resp_signal_short)
plot(t(maxima_loc),resp_signal_short(maxima_loc),'bo');
plot(t(minima_loc),resp_signal_short(minima_loc),'ro');
subplot(313), hold on
plot(phase_bin(resp_signal_steps))
%}
end
function [minima_loc, maxima_loc] = correct_double_extrema(minima_loc,maxima_loc,signal)
% [minima, maxima] = correct_double_extrema(minima,maxima)
%
%
% 18-Feb-2020 l.m.gottwald@amsterdamumc.nl

ii = 1;
n_extrema = min(length(minima_loc),length(maxima_loc));
% CHECK WHICH EXTREMA COMES FIRST
if maxima_loc(1,1) > minima_loc(1,1) % minima first, then maxima
    while ii < n_extrema
        % find next maxima
        next = maxima_loc>minima_loc(ii) & maxima_loc<minima_loc(ii+1);
        if sum(next)>1 % remove multiple maxima
            idx_next_remove = find(signal(maxima_loc(next))~=max(signal(maxima_loc(next)))); %remove all not absolute maxima
            idx_next = find(next);
            idx_next_remove = idx_next(idx_next_remove);
            maxima_loc(idx_next_remove)=[];
        elseif sum(next)==0 % remove multiple minima
            if signal(minima_loc(ii)) > signal(minima_loc(ii+1)) %remove all not absolute minima
                idx_next_remove = ii;
            else
                idx_next_remove = ii+1;
            end
            minima_loc(idx_next_remove)=[];               
        end  
        % find next minima
        next = minima_loc>maxima_loc(ii) & minima_loc<maxima_loc(ii+1);
        if sum(next)>1 % remove multiple minima
            idx_next_remove = find(signal(minima_loc(next))~=min(signal(minima_loc(next)))); %remove all not absolute minima
            idx_next = find(next);
            idx_next_remove = idx_next(idx_next_remove);
            minima_loc(idx_next_remove)=[];
        elseif sum(next)==0 % remove multiple maxima
            if signal(maxima_loc(ii)) > signal(maxima_loc(ii+1)) %remove all not absolute maxima
                idx_next_remove = ii;
            else
                idx_next_remove = ii+1;
            end
            maxima_loc(idx_next_remove)=[];              
        end    
        n_extrema = min(length(minima_loc),length(maxima_loc));
        ii = ii+1;         
    end    
else % maxima first, then minima
    while ii < n_extrema    
        % find next minima
        next = minima_loc>maxima_loc(ii) & minima_loc<maxima_loc(ii+1);
        if sum(next)>1 % remove multiple minima
            idx_next_remove = find(signal(minima_loc(next))~=min(signal(minima_loc(next)))); %remove all not absolute minima
            idx_next = find(next);
            idx_next_remove = idx_next(idx_next_remove);
            minima_loc(idx_next_remove)=[];
        elseif sum(next)==0 % remove multiple maxima
            if signal(maxima_loc(ii)) > signal(maxima_loc(ii+1)) %remove all not absolute maxima
                idx_next_remove = ii;
            else
                idx_next_remove = ii+1;
            end
            maxima_loc(idx_next_remove)=[];            
        end
        
        % find next maxima
        next = maxima_loc>minima_loc(ii) & maxima_loc<minima_loc(ii+1);
        if sum(next)>1 % remove multiple maxima
            idx_next_remove = find(signal(maxima_loc(next))~=max(signal(maxima_loc(next)))); %remove all not absolute maxima
            idx_next = find(next);
            idx_next_remove = idx_next(idx_next_remove);
            maxima_loc(idx_next_remove)=[];
        elseif sum(next)==0 % remove multiple minima
            if signal(minima_loc(ii)) > signal(minima_loc(ii+1)) %remove all not absolute minima
                idx_next_remove = ii;
            else
                idx_next_remove = ii+1;
            end
            minima_loc(idx_next_remove)=[];              
        end        
        n_extrema = min(length(minima_loc),length(maxima_loc));
        ii = ii+1;        
    end
end
end
function phase_bin = phase_binning(signal,maxima, minima, nbins)
% phase_bin = phase_binning(signal,maxima, minima)
% 
% This functions bins SIGNAL with MAXIMA and MINIMA into NBINS
% MAXIMA and MINIMA are locations in SIGNAL
% NBINS are between consecutive MAXIMA and MINIMA
% e.g. nbins=100; 1 to 100, 1 for expiration, 100 for inspiration state
% 
%
% 30-Oct-2019 Lukas Gottwald, l.m.gottwald@amsterdamumc.nl

%% PERFORM PHASE BINNING 
phase_bin = zeros(size(signal));
nmax = length(maxima);
nmin = length(minima);

% CHECK WHICH EXTREMA COMES FIRST
if maxima(1,1) > minima(1,1) % minima first, then maxima
    for ii = 1:min(nmax,nmin)-1
        distance = (maxima(ii,1)-minima(ii,1))+1; % min_n to max_n
        tmp = fix(linspace(1,nbins,distance));
        phase_bin(minima(ii,1):maxima(ii,1))=tmp';
        
        distance = (minima(ii+1,1)-maxima(ii,1))+1; % max_n to min_n+1
        tmp = fliplr(fix(linspace(1,nbins,distance)));
        phase_bin(maxima(ii,1):minima(ii+1,1))=tmp';
    end
    % create smooth boundaries
    phase_bin(1:minima(1,1)-1) = flipud(phase_bin(minima(1,1)+1:minima(1,1)+minima(1,1)-1));
    phase_bin(minima(ii+1,1)+1:end) = flipud(phase_bin(minima(ii+1,1)-(length(phase_bin)-minima(ii+1,1)):minima(ii+1,1)-1));
else % maxima first, then minima
    for ii = 1:min(nmax,nmin)-1
        distance = (minima(ii,1)-maxima(ii,1))+1; % max_n to min_n
        tmp = fliplr(fix(linspace(1,nbins,distance)));
        phase_bin(maxima(ii,1):minima(ii,1))=tmp';
        
        distance = (maxima(ii+1,1)-minima(ii,1))+1; % min_n to max_n+1
        tmp = fix(linspace(1,nbins,distance));
        phase_bin(minima(ii,1):maxima(ii+1,1))=tmp';
    end
    % create smooth boundaries
    phase_bin(1:maxima(1,1)-1) = flipud(phase_bin(maxima(1,1)+1:maxima(1,1)+maxima(1,1)-1));
    phase_bin(maxima(ii+1,1)+1:end) = flipud(phase_bin(maxima(ii+1,1)-(length(phase_bin)-maxima(ii+1,1)):maxima(ii+1,1)-1));    
end
end
   