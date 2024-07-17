function [extr2, nPhases] = retro_vEye_Sort(labels,vEye,nPhases)

%% find vEye in MRecon.Data
typ = labels.Index.typ;

% set typ == 2 to 1
typ(typ==2) = 1;

noisemeas = sum(labels.Index.typ == 5);
phasecorrmeas = sum(labels.Index.typ == 3);
nchan = size(labels.CoilNrsPerStack{1,1},1);

% hack here, set the first readout (1:nchan) of points equal to
% point=nchan+1
vEye(1:nchan) = vEye(nchan+1);

% normalize 0 to 1
vEye = vEye - min(vEye); vEye = vEye/max(vEye(:));

%% Create range of breathing locations, sort into nPhases resp phases
vEye = smooth(vEye,100);
limit_low = min(vEye);%median(vEye) - std(vEye);

% so we want resp phases with roughly:
readoutsPerResp = numel(find(vEye >= limit_low))/nPhases;

step = 0.002; %abs(mean(diff(respRange))/99);

% start with max_beamLoc and add readouts until we get to readoutsPerResp
% the last phase will have the least amount of readouts
extr2 = zeros(size(vEye))+100;               
limits = zeros(2,nPhases);                  % the limits for plotting
for ms = 1:nPhases
    temp_ind = [];
    if ms == 1
        startbeamLoc = max(vEye);
    end
    limits(1,ms) = startbeamLoc;
    if ms < nPhases
        while numel(temp_ind) < readoutsPerResp
            temp_ind = cat(1,temp_ind,find(vEye <= startbeamLoc & ...
                vEye > (startbeamLoc - step)));
            startbeamLoc = startbeamLoc - step;
        end
    else    % in the last resp phase simply find points that haven't been sorted yet
        temp_ind = find(extr2 == 100 & vEye >= limit_low);
        startbeamLoc = limit_low;
    end
    extr2(temp_ind) = ms-1;
    limits(2,ms) = startbeamLoc;
end

%% optional limit plots

c = lines(nPhases);
if 1
    figure(1);  clf;
    subplot 211;
    plot(vEye,'Color','k'); hold on;  
    x = [1 length(vEye)];
    for i = 1:nPhases
        I = patch([x fliplr(x)],[limits(1,i) limits(1,i) fliplr([limits(2,i) limits(2,i)])],'k');
        I.FaceColor = c(i,:); I.FaceAlpha = 0.2;
    end
    
    xlim([1 length(vEye)])
    xlabel('readout number'); ylabel('vEye signal (a.u.)')
    
    subplot 212; hold on
    histogram(vEye,200, 'FaceColor','k');
    for i = 1:nPhases
        area([limits(2,i), limits(1,i)], [max(ylim), max(ylim)],'FaceColor',...
            c(i,:),'FaceAlpha',0.2);
    end
    xlabel('vEye signal (a.u.)'); ylabel('bin count')
end

set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','FontName'),'FontName','Microsoft Yahei')
set(gcf,'Position',[964 53 954 913])