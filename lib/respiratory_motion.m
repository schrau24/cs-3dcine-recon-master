function [Res_Signal, sg_signal_final] = respiratory_motion(mrecon)

% LOAD DATA
SG = mrecon.Copy;
SG.Parameter.Parameter2Read.typ = 1;
SG.Parameter.Labels.Index.typ(SG.Parameter.Labels.Index.typ==2)=1; % include rejected data
SG.Parameter.Parameter2Read.Update;
SG.ReadData;

SG.RandomPhaseCorrection;
SG.RemoveOversampling;
SG.PDACorrection;
SG.DcOffsetCorrection;
SG.MeasPhaseCorrection;

noisemeas = sum(SG.Parameter.Labels.Index.typ == 5);
phasecorrmeas = sum(SG.Parameter.Labels.Index.typ == 3);
data_in_use_index = false(size(SG.Parameter.Labels.Index.rtop));
data_in_use_index(1+noisemeas+phasecorrmeas:SG.Parameter.Labels.OriginalLabelLength)=true;

SG.Data = abs(SG.Data);

%% Define SG points to be used
% channels
chan = SG.Parameter.Labels.CoilNrs(:,1);
tmp = SG.Parameter.Labels.Index.ky(data_in_use_index)==0 & SG.Parameter.Labels.Index.kz(data_in_use_index) == 0 & SG.Parameter.Labels.Index.chan(data_in_use_index) == chan(1);
sg_sampling_index = find(tmp); clear tmp;

% sort for channels
clear signal
for ii=1:size(chan,1)
    signal(:,:,ii) = SG.Data(:,SG.Parameter.Labels.Index.ky(data_in_use_index)==0 & SG.Parameter.Labels.Index.kz(data_in_use_index) == 0 & SG.Parameter.Labels.Index.chan(data_in_use_index) == chan(ii));
end
% The second dimension is now the time dimension for breathing

% Respiratory motion detection
ZIP = signal;
%ZIP = abs(ZIP);
ZIP = abs(fftshift(ifft(ZIP,400,1),1));
n1 = size(ZIP,2);

%Normalization of each projection in each coil element
ZIP=ProjNorm(ZIP);%Normalization includes temporal smoothing
figure,imagesc(abs(ZIP(:,:,10))),axis image, axis off, colormap(gray),title('Respiratory Motion')


% STEP 1: find the coil elements with good representation of respiratory motion
%         from the late enhancement spokes
[Coil,Res_Signal_Post] = MC_Step1(ZIP,n1);

% STEP 2: Estimate motion signal using PCA from the concatated coil elements
%         Those coil elements were selected in the first step
[SI,corrm,Res_Signal,ZIP1] = MC_Step2(ZIP,Coil,n1,Res_Signal_Post);

% check to flip resp signal if expiration is on bottom
if sum(Res_Signal < 0.5) > sum(Res_Signal > 0.5)
    Res_Signal = -Res_Signal + 1;
end

figure(1);clf;
plot(Res_Signal(:)*100+220,'r')
axis image
xlabel('view number')
ylabel('resp signal (a.u.)')
title('Respiratory Motion')

%% WRITE LABEL
sg_signal = zeros(size(data_in_use_index));
for ii = 1:size(Res_Signal,1)
    if ii == 1
        sg_signal(1:sg_sampling_index(ii+1)-1) = Res_Signal(ii,1);
    elseif ii < size(Res_Signal,1)
        sg_signal(sg_sampling_index(ii):sg_sampling_index(ii+1)-1) = Res_Signal(ii,1);
    else
        sg_signal(sg_sampling_index(ii):end) = Res_Signal(ii,1);
    end
end

sg_signal_final = sg_signal;
