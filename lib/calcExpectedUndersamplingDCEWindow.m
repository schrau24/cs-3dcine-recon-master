clc;
TR = 3.7;%RESP.Parameter.Labels.RepetitionTime;  % in ms

% using given dceWindow for DCE, this would give:
dceWindow = 3; 
nReadoutsPerWindow = round(dceWindow*1000/TR);

% PROUDSampling = load('/home/emschrauben/lood_storage/divi/Projects/dcerecon/Spirals_for_DCE/01-06-2021/001_010621_Y192_Z46_scanTime5.00_mode1.0-0.0_NS48422_Np100_t1_flip1.dat');
[fName, pName] = uigetfile('*.dat','select PROUD sampling file');
PROUDSampling = load(fullfile(pName,fName));
sy = PROUDSampling(:,1); sy = sy-min(sy)+1;
sz = PROUDSampling(:,2); sz = sz-min(sz)+1;
lsy = length(min(sy):max(sy));
lsz = length(min(sz):max(sz));


% find those unique locations and set first nReadoutsPerWindow = 1
unq = unique(cat(2,sy(1:nReadoutsPerWindow), sz(1:nReadoutsPerWindow)),'rows');
idx = sub2ind([lsy,lsz],unq(:,1),unq(:,2));
tempMask = zeros(lsy,lsz);
tempMask(idx) = 1;

mask = repmat(permute(tempMask,[3 1 2]),[lsy, 1, 1]);

R = numel(mask)/sum(mask(:)) /4*pi;

disp(['DCE window = ' num2str(dceWindow) ' s and TR = ' num2str(round(TR,2))...
    ' ms --> R = ' num2str(round(R,2))])
figure(1); clf; imshow(tempMask)