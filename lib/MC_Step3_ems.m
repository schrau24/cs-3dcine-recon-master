function [Res_Signal2]=MC_Step3_ems(Res_Signal)
% Motion detection: step 3

% find the pre- and post-contrast peaks, subtract the mean of peaks for each
[~, Peak_Index]=findpeaks(double(Res_Signal),'MinPeakProminence',0.02);
% Peak_Index=find(islocalmin(double(Res_Signal),'MinProminence',0.02));
figure,plot(Res_Signal);
hold on;
plot(Peak_Index,Res_Signal(Peak_Index),'ro');

if length(Peak_Index) > 1
    
    % find peak diff outliers, the injection phase will be between 100 and 350
    pkDiffsOutliers = Peak_Index(isoutlier(abs(diff(Res_Signal(Peak_Index)))));
    if length(pkDiffsOutliers) == 1 % add the next peak as well
        pkDiffsOutliers(2) = Peak_Index(find(Peak_Index>pkDiffsOutliers,1,'first'));
    end
    pkDiffsOutliers(pkDiffsOutliers<100 | pkDiffsOutliers > 350) = [];
    
    preRes = [zeros(10,1); Res_Signal(10:pkDiffsOutliers(1)-1)];  % again skip first 10 points
    preRes = preRes - mean(Res_Signal(Peak_Index(1):pkDiffsOutliers(1)));
    
    postRes = Res_Signal(pkDiffsOutliers(end):end);
    p = polyfit(1:length(postRes),postRes,2);
    postRes = postRes - polyval(p,1:length(postRes))';
    
    % during injection, fit
    midInjection = Res_Signal(pkDiffsOutliers(1)+1:pkDiffsOutliers(end)-1);
    p = polyfit(1:length(midInjection),midInjection,1);
    midInjection = midInjection - polyval(p,1:length(midInjection))';
    
    Res_Signal2 = cat(1,preRes,midInjection,postRes) + 1;
    
end

figure; plot(Res_Signal2)