function [extr2, nPhases] = retro_SG_RegistrationSort(labels,sg,nPhases)

% extr2 = retro_resp_binning(labels,sg)
%
% retrospective respiratory binning (for one respiratory frame)
% labels: labels
% vEye: -double(mrecon.Parameter.Labels.Index.na)
%
% 26-Nov-2020 e.m.schrauben@amsterdamumc.nl

%% find sg in MRecon.Data
typ = labels.Index.typ;

% set typ == 2 to 1
typ(typ==2) = 1;

noisemeas = sum(labels.Index.typ == 5);
phasecorrmeas = sum(labels.Index.typ == 3);
nchan = size(labels.CoilNrsPerStack{1,1},1);

% hack here, set the first readout (1:nchan) of points equal to
% point=nchan+1
sg(1:nchan) = sg(nchan+1);

% normalize 0 to 1
sg = sg - min(sg); sg = sg/max(sg(:));

%% Create range of breathing locations, sort into nPhases resp phases

% here we take the the maximum beam_loc as the high limit, and create a low
% limit based on the median - std
limit_low = min(sg);%median(sg) - std(sg);
% 
% % manually set nPhases = 6
% nPhases = 6;

% so we want resp phases with roughly:
readoutsPerResp = numel(find(sg >= limit_low))/nPhases;

step = 0.002; %abs(mean(diff(respRange))/99);

% start with max_beamLoc and add readouts until we get to readoutsPerResp
% the last phase will have the least amount of readouts
extr2 = zeros(size(sg))+100;               % initiate all at 100 for MRecon parameter update later
limits = zeros(2,nPhases);                  % the limits for plotting
for ms = 1:nPhases
    temp_ind = [];
    if ms == 1
        startbeamLoc = max(sg);
    end
    limits(1,ms) = startbeamLoc;
    if ms < nPhases
        while numel(temp_ind) < readoutsPerResp
            temp_ind = cat(1,temp_ind,find(sg <= startbeamLoc & sg > (startbeamLoc - step)));
            startbeamLoc = startbeamLoc - step;
            if startbeamLoc < min(sg)
                sprintf('Sorting of sg signal is not going well')
                break
            end
        end
    else    % in the last resp phase simply find points that haven't been sorted yet
        temp_ind = find(extr2 == 100 & sg >= limit_low);
        startbeamLoc = limit_low;
    end
    extr2(temp_ind) = ms-1;
    limits(2,ms) = startbeamLoc;
end

%% optional limit plots

c = lines(nPhases);
if 1
    figure(2);  clf;
    subplot 211;
    plot(sg,'Color','k'); hold on;  
    x = [1 length(sg)];
    for i = 1:nPhases
        I = patch([x fliplr(x)],[limits(1,i) limits(1,i) fliplr([limits(2,i) limits(2,i)])],'k');
        I.FaceColor = c(i,:); I.FaceAlpha = 0.2;
    end
    
    xlim([1 length(sg)])
    xlabel('readout number'); ylabel('sg signal (a.u.)')
    
    subplot 212; hold on
    histogram(sg,200, 'FaceColor','k');
    for i = 1:nPhases
        area([limits(2,i), limits(1,i)], [max(ylim), max(ylim)],'FaceColor',...
            c(i,:),'FaceAlpha',0.2);
    end
    xlabel('sg signal (a.u.)'); ylabel('bin count')
end

set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','FontName'),'FontName','Microsoft Yahei')
set(gcf,'Position',[964 53 954 913])