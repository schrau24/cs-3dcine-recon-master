function [extr2, nPhases] = retro_VitalEye_RegistrationSort(mrecon, ...
                    nPhases);
                
 SG = mrecon.Copy;
 labels = SG.Parameter.Labels;
 VEye = SG.Parameter.Labels.Index.na;
 typ = labels.Index.typ;
 
 
noisemeas = sum(labels.Index.typ == 5);
phasecorrmeas = sum(labels.Index.typ == 3);
nchan = size(labels.CoilNrsPerStack{1,1},1);

 %Select only the readouts not the noise meas and phasecorr, add at the end
 VEye = VEye(typ==1);
 VEye = int64(VEye);
 VEye(VEye > 20000) = VEye(VEye > 20000) -max(VEye);
 VEye = double(VEye);
 
 VEye = -1*VEye;
 VEye = VEye - min(VEye); VEye = VEye/max(VEye(:));
 
 %% Create range of breathing locations, sort into nPhases resp phases

% here we take the the maximum beam_loc as the high limit, and create a low
% limit based on the median - std
limit_low = min(VEye);%median(VEye) - std(VEye);

% so we want resp phases with roughly:
readoutsPerResp = numel(find(VEye >= limit_low))/nPhases;

step = 0.005; %abs(mean(diff(respRange))/100);

% start with max_beamLoc and add readouts until we get to readoutsPerResp
% the last phase will have the least amount of readouts
extr2 = zeros(size(VEye))+100;               % initiate all at 100 for MRecon parameter update later
limits = zeros(2,nPhases);                  % the limits for plotting
for ms = 1:nPhases
    temp_ind = [];
    if ms == 1
        startbeamLoc = max(VEye);
    end
    limits(1,ms) = startbeamLoc;
    if ms < nPhases
        while numel(temp_ind) < readoutsPerResp
            temp_ind = cat(1,temp_ind,find(VEye <= startbeamLoc & ...
                VEye > (startbeamLoc - step)));
            startbeamLoc = startbeamLoc - step;
        end
    else    % in the last resp phase simply find points that haven't been sorted yet
        temp_ind = find(extr2 == 100 & VEye >= limit_low);
        startbeamLoc = limit_low;
    end
    extr2(temp_ind) = ms-1;
    limits(2,ms) = startbeamLoc;
end

% finally add the noise and phase corr points back in
extr2 = cat(1,100*ones(noisemeas+phasecorrmeas,1), extr2);
%% optional limit plots

c = lines(nPhases);
if 1
    figure(1);  clf;
    subplot 211;
    plot(VEye,'Color','k'); hold on;  
    x = [1 length(VEye)];
    for i = 1:nPhases
        I = patch([x fliplr(x)],[limits(1,i) limits(1,i) fliplr([limits(2,i) limits(2,i)])],'k');
        I.FaceColor = c(i,:); I.FaceAlpha = 0.2;
    end
    
    xlim([1 length(VEye)])
    xlabel('readout number'); ylabel('VEye signal (a.u.)')
    
    subplot 212; hold on
    histogram(VEye,200, 'FaceColor','k');
    for i = 1:nPhases
        area([limits(2,i), limits(1,i)], [max(ylim), max(ylim)],'FaceColor',...
            c(i,:),'FaceAlpha',0.2);
    end
    xlabel('VEye signal (a.u.)'); ylabel('bin count')
end

set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','FontName'),'FontName','Microsoft Yahei')
set(gcf,'Position',[964 53 954 913])
end