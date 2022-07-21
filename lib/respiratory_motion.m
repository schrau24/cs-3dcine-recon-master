function [Res_Signal, sg_signal_final] = respiratory_motion(mrecon, t1Flag)

% LOAD DATA
SG = mrecon.Copy;
SG.Parameter.Parameter2Read.typ = 1;
SG.Parameter.Labels.Index.typ(SG.Parameter.Labels.Index.typ==2)=1; % include rejected data
SG.Parameter.Parameter2Read.Update;
SG.ReadData;
% SG.K2IM; % Fourier transform in frequency encoding direction 

noisemeas = sum(SG.Parameter.Labels.Index.typ == 5);
phasecorrmeas = sum(SG.Parameter.Labels.Index.typ == 3);
data_in_use_index = false(size(SG.Parameter.Labels.Index.rtop));
data_in_use_index(1+noisemeas+phasecorrmeas:SG.Parameter.Labels.OriginalLabelLength)=true;  %rm: you already only read in data of typ ==1?

SG.Data = abs(SG.Data);

%% Define SG points to be used
% channels
chan=SG.Parameter.Labels.CoilNrs(:,1);
% flow encoding directions
extr1=SG.Parameter.Parameter2Read.extr1;
tmp =SG.Parameter.Labels.Index.ky(data_in_use_index)==0 & SG.Parameter.Labels.Index.kz(data_in_use_index) == 0 & SG.Parameter.Labels.Index.chan(data_in_use_index) == chan(1) & SG.Parameter.Labels.Index.extr1(data_in_use_index)== extr1(1);
sg_sampling_index = find(tmp);

% sort for channels and flow encodings
clear signal
for ii=1:size(chan,1)
    for kk=1:size(extr1,1)
        signal(:,:,ii,kk)=SG.Data(:,SG.Parameter.Labels.Index.ky(data_in_use_index)==0 & SG.Parameter.Labels.Index.kz(data_in_use_index) == 0 & SG.Parameter.Labels.Index.chan(data_in_use_index) == chan(ii) & SG.Parameter.Labels.Index.extr1(data_in_use_index)== extr1(kk)); % use if mean is used above
    end
end
% The secon dimension is now the time dimension for breathing

% Respiratory motion detection
ZIP = signal;
ZIP = abs(ZIP);
%Normalization of each projection in each coil element

% the last 2/3 of the spokes were used as the late enhancement phase for motion detection
% as described in the paper, if t1, use all data
n1 = round(size(ZIP,2)*2/3);
if t1Flag
    n1 = size(ZIP,2);
end
% STEP 1: find the coil elements with good representation of respiratory motion
%         from the late enhancement spokes
[Coil,Res_Signal_Post] = MC_Step1(ZIP,n1);

%STEP 2: Estimate motion signal using PCA from the concatated coil elements
%Those coil elements were selected in the first step
[SI,corrm,Res_Signal,ZIP1] = MC_Step2(ZIP,Coil,n1,Res_Signal_Post);

% Eric, skip this step if no injection or if t1 mapping
% %Step 3: You noticed that the signal is not flat, due to the contrast
% %injection. So, now let's estimate the envelop of the signal and substract it
if ~t1Flag
    Res_Signal = MC_Step3_ems(Res_Signal);
end
close all

% find beginning outliers and set to mean to avoid strange sorting
outLiers = isoutlier(Res_Signal);
if length(tmp) > 1
    Res_Signal(outLiers(1:10)) = mean(Res_Signal(11:end));
end
% check to flip resp signal if expiration is on bottom
if sum(Res_Signal < 0.5) > sum(Res_Signal > 0.5)
    Res_Signal = -Res_Signal + 1;
end

figure(1);clf;
plot(Res_Signal(:)*100+220,'r')
xlabel('view number')
ylabel('resp signal (a.u.)')
title('Respiratory Motion')

%% WRITE LABEL
% sg_signal = zeros(size(data_in_use_index));
sg_signal = zeros(SG.Parameter.Labels.OriginalLabelLength,1);
for ii=1:size(Res_Signal,1)
    if ii==1
        sg_signal(1:sg_sampling_index(ii+1)-1)=Res_Signal(ii,1);
    elseif ii < size(Res_Signal,1)
        sg_signal(sg_sampling_index(ii):sg_sampling_index(ii+1)-1)=Res_Signal(ii,1);
    else
        sg_signal(sg_sampling_index(ii):end)=Res_Signal(ii,1);
    end
end


sg_signal_final = sg_signal;
