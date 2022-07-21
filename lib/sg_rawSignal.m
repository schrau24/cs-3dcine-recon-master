function sg_signal_final = sg_rawSignal(mrecon)%% SELF GATIG TEST SCRIPT


% LOAD DATA
SG = mrecon.Copy;
SG.Parameter.Parameter2Read.typ = 1;
SG.Parameter.Labels.Index.typ(SG.Parameter.Labels.Index.typ==2)=1; % include rejected data
SG.Parameter.Parameter2Read.Update;
SG.ReadData;
SG.K2IM; % Fourier transform in frequency encoding direction 

noisemeas = sum(SG.Parameter.Labels.Index.typ == 5);
phasecorrmeas = sum(SG.Parameter.Labels.Index.typ == 3);
data_in_use_index = false(size(SG.Parameter.Labels.Index.rtop));
data_in_use_index(1+noisemeas+phasecorrmeas:SG.Parameter.Labels.OriginalLabelLength)=true;  %rm: you already only read in data of typ ==1?

SG.Data = abs(SG.Data);

%% Sort data according to coils and flow encodings (extr1 dimension)
% channels
chan=SG.Parameter.Labels.CoilNrs(:,1);

% flow encoding directions
extr1=SG.Parameter.Parameter2Read.extr1;

% sort for channels and flow encodings
clear signal
for ii=1:size(chan,1)
    for kk=1:size(extr1,1)
        signal(:,:,ii,kk)=SG.Data(:,SG.Parameter.Labels.Index.ky(data_in_use_index)==0 & SG.Parameter.Labels.Index.kz(data_in_use_index) == 0 & SG.Parameter.Labels.Index.chan(data_in_use_index) == chan(ii) & SG.Parameter.Labels.Index.extr1(data_in_use_index)== extr1(kk)); % use if mean is used above
    end
end
% The secon dimension is now the time dimension for breathing

%% SVD analysis for all coils
clear resp_coils
signal_center_coil=mean(signal,4); % average over flow encoding

% perform singular value decomposition for (readout x time) for all coils - uh so many coils!
for i=1:size(chan,1)
    [U,S,V]=svd(signal_center_coil(:,:,i)'); %take transpose of signal_center_coil so time series is stored in columns of U
    resp_coils(:,i)=U(:,1);  % better use no smoothing here
end

%% Define SG points to be used

tmp =SG.Parameter.Labels.Index.ky(data_in_use_index)==0 & SG.Parameter.Labels.Index.kz(data_in_use_index) == 0 & SG.Parameter.Labels.Index.chan(data_in_use_index) == chan(1) & SG.Parameter.Labels.Index.extr1(data_in_use_index)== extr1(1);
sg_sampling_index = find(tmp);
sg_sampling_diff=diff(sg_sampling_index);

%% Perform bandpass filter
sampling_t = mean(sg_sampling_diff)/length(chan)*SG.Parameter.Scan.TR/1000; %[s]
sampling_freq = 1/sampling_t; % [Hz]
f1= 3/60; % 10 respirations per minute [Hz]
f2= 40/60; % 40 respirations per minute [Hz]
input = resp_coils;
clear sg_signal_bp
for ii=1:size(chan,1)
    sg_signal_bp(:,ii) = bpfilt(input(:,ii)',f1,f2,sampling_freq,0);
end

figure(43); clf; %for the method comparing coils we also plot all coils in figure 40
c = jet(size(signal,3));
for j=1:size(signal,3) 
    plot(sg_signal_bp(:,j)+j/100000,'Color',c(j,:)); hold on
end
title('plot of signal of all coils'); 
set(figure(43),'Units','Normalized');
set(figure(43),'Position',[0    0.0454    1.0000    0.8481])
legend(string(1:size(signal,3)))

%% Define coils to be used and average signal
%% generate prompt
promptcell = {...
    'Main coil:'....
    };
dlg_title = 'Pick Main coil';
num_lines = 1;
defaultans = {...
    num2str(11)
    };
answer          = inputdlg(promptcell,dlg_title,num_lines,defaultans);
maincoil   = str2double(answer{1,1});


% find coils that correlate the most (correlation and anti-correlation of coils)
th = 0.9; % threshold
R=corrcoef(sg_signal_bp);
R(abs(R)<th)=0;
R(R==1)=0;
% figure,imagesc(R)
corrcoils_neg = R(:,maincoil)<-th;
corrcoils_pos = R(:,maincoil)>th; 
corrcoils_pos(maincoil) = 1;  %renske added, you also want to consider the coil with highest signal itself in further analysis.  
chan_used = find(corrcoils_pos | corrcoils_neg);
% smooth over period of < 4 se

sg_signal_bp_mean = smooth(mean([sg_signal_bp(:,corrcoils_pos),-sg_signal_bp(:,corrcoils_neg)],2),sampling_freq); % smooth signal with neighbors
% figure, plot(sg_signal_bp_mean)

% 
% sg_signal_bp_mean = sg_signal_bp(:,[3]);
% sg_signal_bp_mean = smooth(sg_signal_bp_mean,floor(4/sampling_t)); % SMA over <3 sec

%% WRITE LABEL
sg_signal = zeros(size(data_in_use_index));
for ii=1:size(sg_signal_bp_mean,1)
    if ii==1
        sg_signal(1:sg_sampling_index(ii+1)-1)=sg_signal_bp_mean(ii,1);
    elseif ii < size(sg_signal_bp_mean,1)
        sg_signal(sg_sampling_index(ii):sg_sampling_index(ii+1)-1)=sg_signal_bp_mean(ii,1);
    else
        sg_signal(sg_sampling_index(ii):end)=sg_signal_bp_mean(ii,1);
    end
end


sg_signal_final = sg_signal;

