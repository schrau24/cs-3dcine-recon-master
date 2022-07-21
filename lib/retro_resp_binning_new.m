function labels = retro_resp_binning_new(labels,method)
% labels = retro_resp_binning_new(labels)
%
% retrospective respiratory binning
% labels: mrecon.Parameter.Labels
% method: method == 0   RNAV binning
% method: method == 1   VitalEye binning
% method: method == 2   self gating
%
% 10-Jul-2019 l.m.gottwald@amsterdamumc.nl

%% get respiratory data
nchan = size(labels.CoilNrsPerStack{1,1},1);
nenc = max(labels.Index.extr1)+1;
noisemeas = sum(labels.Index.typ == 5);
phasecorrmeas = sum(labels.Index.typ == 3);
data_in_use_index(1+noisemeas+phasecorrmeas:labels.OriginalLabelLength)=true;

if method == 0 % RNAV
    resp_signal = double(labels.Index.rnav(data_in_use_index));
elseif method == 1 % VitalEye
    resp_signal = -double(labels.Index.na(data_in_use_index));
elseif method == 2 % SG
    resp_signal = double(labels.Index.sg(data_in_use_index));    
else
    error('ERROR in retro_resp_binning: no method defined.');
end

%% Normalise values of an array to be between -1 and 1
% delete leading zeros
startindex = find(resp_signal,1,'first');
resp_signal(1:startindex)=resp_signal(startindex);

max_range_value = max(resp_signal(resp_signal~=0));
min_range_value = min(resp_signal(resp_signal~=0));
resp_signal_norm = zeros(size(resp_signal));
resp_signal_norm(resp_signal~=0) = -1 + 2 .* (resp_signal(resp_signal~=0)-min_range_value) ./ (max_range_value - min_range_value);

% resp_signal_norm same for all channel per stack: pick most frequent value in a sample
step = double(nchan*nenc);
resp_signal_point = [1:step:size(resp_signal_norm,1)]';
resp_signal_norm = resp_signal_norm(1:step:end);

% smooth signal
resp_signal_norm = smooth(resp_signal_norm,100); %SMA100

%% Find peaks
mpp = mean(abs(resp_signal_norm))/4;
mpd = 100;

[max_values,max_locs, ~, ~] = findpeaks(resp_signal_norm, 'MinPeakProminence', mpp,'MinPeakDistance',mpd);
[min_values,min_locs, ~, ~] = findpeaks(-resp_signal_norm, 'MinPeakProminence', mpp,'MinPeakDistance',mpd);

max_values_mean = mean(max_values,'omitnan');
min_values_mean = mean(-min_values,'omitnan');
max_values_sd = std(max_values,'omitnan');
min_values_sd = std(min_values,'omitnan');
minmax_range = max_values_mean-min_values_mean;

%% PHASE BINNING
% maxima = 1 = expiration
% minima = 0 = inspiration
nbins = 100;
phase_bin = phase_binning(resp_signal_norm,max_locs, min_locs, nbins);  

resp_signal_typ = false(size(resp_signal_norm)); % by default rejected
resp_signal_bin = zeros(size(resp_signal_norm)); % by default bin #0

limits=[];
limits(1)=100; % expiration state
limits(2)=40; % 60% acceptance rate / rejected
% limits(3)=30; % 70% acceptance rate / rejected

for ilimit = 1:size(limits,2)-1
    idx = phase_bin<=limits(ilimit) & phase_bin>limits(ilimit+1);
    resp_signal_typ(idx)=true;
%     resp_signal_bin(idx)= ilimit;
end
% % TESTPLOT
% figure,hold on
% title('RESP SIGNAL');
% plot(1:length(resp_signal_norm),resp_signal_norm);
% plot(max_locs,resp_signal_norm(max_locs),'go');
% plot(min_locs,resp_signal_norm(min_locs),'ro');

% plot(1:length(resp_signal_bin),resp_signal_bin)
% 
% subplot(212), hold on
% title('VE');
% plot(1:length(ve),ve);
% plot(maxlocs_ve,ve(maxlocs_ve),'go');
% plot(minlocs_ve,ve(minlocs_ve),'ro');



%% AMPLITUDE BINNING MIN95 METHOD
% % Define acceptance window(s) and bin respiratory signal
% limits=[];
% % % mean peak + 2SD range % LG: these limits should be defined smarter in the future
% % limits(1) = max_values_mean+2*max_values_sd;
% % limits(2) = max_values_mean+2*max_values_sd-minmax_range*0.5;
% % limits(3) = max_values_mean+2*max_values_sd-minmax_range*0.8;
% 
% % MIN95 method
% limits(1) = quantile(resp_signal_norm,0.975);
% limits(2) = quantile(resp_signal_norm,0.575);
% 
% % % with R=8 maximal acceleration
% % rmax = 8;
% % labels_full = double(labels.YResolution) * double(labels.ZResolution)*double(labels.KzOversampleFactor) * double(range(labels.Index.card)+1) * double(range(labels.Index.extr1)+1) * double(length(labels.CoilNrs));
% % labels_rmax = fix(labels_full / rmax);
% % resp_signal_norm_sorted = sort(resp_signal_norm(resp_signal_norm<limits(1)));
% % limits(2) = resp_signal_norm_sorted(length(resp_signal_norm_sorted)-labels_rmax);
% 
% % limits(end+1) = quantile(resp_signal_norm,0.875)
% 
% resp_signal_typ = false(size(resp_signal_norm)); % by default rejected
% resp_signal_bin = zeros(size(resp_signal_norm)); % by default bin #0
% 
% for ii = 1:size(limits,2)-1
%     idx = resp_signal_norm<limits(ii) & resp_signal_norm>limits(ii+1);
%     resp_signal_typ(idx)=true;
%     resp_signal_bin(idx)= ii-1;
% end
%% Rescale
resp_signal_norm_new=resp_signal;
phase_bin_new=resp_signal;
resp_signal_bin_new=resp_signal;
resp_signal_typ_new=resp_signal;
resp_signal_point(length(resp_signal_point)+1)=length(resp_signal);
for ii=1:size(resp_signal_point,1)-1
    resp_signal_norm_new(resp_signal_point(ii):resp_signal_point(ii+1)-1)=resp_signal_norm((ii));
    phase_bin_new(resp_signal_point(ii):resp_signal_point(ii+1)-1)=phase_bin((ii));
    resp_signal_bin_new(resp_signal_point(ii):resp_signal_point(ii+1)-1)=resp_signal_bin((ii));
    resp_signal_typ_new(resp_signal_point(ii):resp_signal_point(ii+1)-1)=resp_signal_typ((ii));
end
resp_signal_norm = resp_signal_norm_new;
phase_bin = phase_bin_new;
resp_signal_bin = resp_signal_bin_new;
resp_signal_typ = resp_signal_typ_new;
clear resp_signal_norm_new phase_bin_new resp_signal_bin_new resp_signal_typ_new

%% Overwrite labels
typ = labels.Index.typ(data_in_use_index);
typ(~resp_signal_typ) = 2; % label data out of respiratory limits as rejected
labels.Index.typ(data_in_use_index) = typ;
labels.Index.extr2(data_in_use_index) = uint16(resp_signal_bin); % LG: use extr2 dimension for resp binning

labels.Index.resp_phase = zeros(size(labels.Index.typ));
labels.Index.resp_phase(data_in_use_index) = phase_bin;
labels.Index.resp_signal = zeros(size(labels.Index.typ));
labels.Index.resp_signal(data_in_use_index) = resp_signal_norm;
%% Testplot
% figure, hold on
% plot(1:length(resp_signal_norm),resp_signal_norm,'-k');
% legendlist{1} = 'Respiratory signal';
% plotcolor = winter(size(limits,2));
% plotcolor(end,:) = [1 0 0];
% for ii = 1:size(limits,2)
%     plot(ones(1,length(resp_signal_norm))*(limits(ii)),'-','color',plotcolor(ii,:));
%     if ii == 1
%         legendlist{ii+1} = sprintf('Bin #%d | Upper limit',ii);
%     elseif ii == size(limits,2)
%         legendlist{ii+1} = 'Lower limit';
%     else
%         legendlist{ii+1} = sprintf('Bin #%d',ii);
%     end
% end
% legend(legendlist,'Location','SouthEast');


%     [x,y] = find(resp_signal_typ);
%     plot(x,y*limits(2),'g.');
%     [x,y] = find(~resp_signal_typ);
%     plot(x,y*limits(2),'r.');
% plot(max_locs,resp_signal_norm(max_locs),'go');
% plot(min_locs,resp_signal_norm(min_locs),'ro');
% plot(ones(1,length(resp_signal_norm))*(limits(1)),'k-');
% plot(ones(1,length(resp_signal_norm))*(min_values_mean-1*min_values_sd),'k-');
% plot(ones(1,length(resp_signal_norm))*(limits(2)),'b-');
% xlim([0 100000])
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
