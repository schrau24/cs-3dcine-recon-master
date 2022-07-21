function [Res_Signal2]=MC_Step3_ems_usingSplineFit(Res_Signal)
% Motion detection: step 3
%%
[~, Peak_Index]=findpeaks(double(Res_Signal),'MinPeakProminence',0.02);
figure,plot(Res_Signal);
hold on;
plot(Peak_Index,Res_Signal(Peak_Index),'ro');

if length(Peak_Index) > 1
    %Do a fitting and subtract the fitted signal
    [xData, yData] = prepareCurveData(Peak_Index,Res_Signal(Peak_Index));
    % interp params
    x = 1:length(xData);
    xq = linspace(1,length(xData),length(Res_Signal));
    METHOD = 'spline';
    fitResult = interp1(x,yData,xq,METHOD);
    
    Res_Signal2=Res_Signal-(fitResult'-1);
    
end

% Res_Signal=single(Res_Signal);

figure; plot(Res_Signal2)