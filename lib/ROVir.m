function data = ROVir( data, RP )

% use Regionally Optimized Virtual (ROVir) coil reduction
% DOI: 10.1002/mrm.28706
% this is either a predefined mask stored in RP.data_dir, or a
% spherical/ellipsoidal mask about 1/4 of the FOV in the center of the FOV

% inputs:
% data:         the original sorted k-space data
% RP:           the .pro.dat file with recon info

%% Re-load and perform simple FFT over time-averaged data
MRtemp = MRecon(fullfile(RP.data_dir,RP.data_target));
MRtemp.Parameter.Cardiac.RetroPhases = 1;
MRtemp.Parameter.Parameter2Read.typ = 1;
MRtemp.Parameter.Parameter2Read.Update;
MRtemp.ReadData;

% perform MRecon reconstruction steps
MRtemp.RandomPhaseCorrection;
MRtemp.RemoveOversampling;
MRtemp.PDACorrection;
MRtemp.DcOffsetCorrection;
MRtemp.MeasPhaseCorrection;
MRtemp.Parameter.Parameter2Read.Update;
MRtemp.SortData;
MRtemp.GridData;

% pre-whiten
MRn=MRecon(fullfile(RP.data_dir,RP.data_target));
MRn.Parameter.Recon.ArrayCompression = MRtemp.Parameter.Recon.ArrayCompression;
MRn.Parameter.Recon.ACNrVirtualChannels = MRtemp.Parameter.Recon.ACNrVirtualChannels;
MRn.Parameter.Parameter2Read.typ=5;
MRn.ReadData;
eta=MRn.Data;

Ncoils=size(MRtemp.Data,4);
Nsamples=numel(MRtemp.Data)/Ncoils;

psi = (1/(Nsamples-1))*(eta' * eta);
L = chol(psi,'lower');
L_inv = (inv(L));
L_inv=diag(diag(L_inv)); %using only diagonal values

MRtemp.Data=permute(MRtemp.Data, [1:3 5:length(size(MRtemp.Data)) 4]);
sizeMRDATA = size(MRtemp.Data);
MRtemp.Data=reshape(MRtemp.Data,[Nsamples,Ncoils]);
MRtemp.Data=MRtemp.Data.';
MRtemp.Data = conj(L_inv) * MRtemp.Data;
MRtemp.Data=MRtemp.Data.';
MRtemp.Data=reshape(MRtemp.Data,sizeMRDATA);
MRtemp.Data=ipermute(MRtemp.Data, [1:3 5:length(size(MRtemp.Data)) 4]);
clearvars L_inv L psi MRn eta

MRtemp.RingingFilter;
orig = MRtemp.Copy;
MRtemp.K2I;

img_orig = MRtemp.Data;
MRtemp.CombineCoils;
img_orig_cc = MRtemp.Data/max(MRtemp.Data(:));
clear MRtemp;

%% masking
[x, y, z] = size(img_orig,1:3);

try % try to load signal and interference masks, if not available, use center spherical FOV
    load(fullfile(RP.data_dir,'ROVirMask.mat'))
    
catch   % spherical/ellipsoidal mask
    [X,Y,Z] = ndgrid(1:x,1:y,1:z);
    x = round(x/2); y = round(y/2); z = round(z/2);
    % define ellipsoid regions as 1/4 FOV
    Rz = round(z/4);
    Ry = round(y/4);
    Rx = round(x/4);
    signal = sqrt( ((z-Z)/(Rz)).^2 + ((y-Y)/(Ry)).^2 + ((x-X)/(Rx)).^2 ) <= 1;
    % use a spherical strel to buffer the regions between signal and
    % interference
    interference = ones(size(signal)) - imdilate(signal,strel('sphere',5));
end


%% perform ROVir
% Calculate eigenvector weights
[eigenvec, SIR] = ROVir_calculate_weights(img_orig, logical(signal), logical(interference), true);

% Apply weights to original data
data = permute(data, [4 1:3 5:length(size(data))]);
data_orig = permute(orig.Data, [4 1:3]);
sizeMRDATA = size(data);

% Form virtual coil data (eq. 2)
nc = sizeMRDATA(1);
nv = RP.nROVirCoils;
weights = eigenvec(:, 1:nv);
k_virt = zeros([nv sizeMRDATA(2:end)], 'single');
k_virt_orig = zeros([nv sizeMRDATA(2:4)]);
for ch_j = 1:nv
    disp('... virtual coil '+string(ch_j)+' of '+string(nv))
    for ch_l = 1:nc
        k_virt(ch_j,:,:,:,1,:,1,1,1,1,:) = k_virt(ch_j,:,:,:,1,:,1,1,1,1,:) + ...
            data(ch_l,:,:,:,1,:,1,1,1,1,:) * weights(ch_l,ch_j);
        k_virt_orig(ch_j,:,:,:) = k_virt_orig(ch_j,:,:,:) + data_orig(ch_l,:,:,:) * weights(ch_l,ch_j);
    end
end

data = permute(k_virt, [2:4 1 5:length(size(data))]);  % this is output!

% re-calculate time-averaged recon with ROVir
orig.Data = permute(single(k_virt_orig), [2 3 4 1]);
orig.K2I; orig.CombineCoils;

%% show the original and ROVir time-averaged image

img_rovir_cc = abs(orig.Data)/max(abs(orig.Data(:)));
% find central slice of the signal region for plotting
[~, ~, zz] = ind2sub(size(signal),find(signal));
sl = min(zz) + round(3*range(zz)/4);
origImg = img_orig_cc(:,:,sl);
figure(22); clf;
subplot 131
imshow(origImg,[0 0.2]);
green = cat(3, zeros(size(origImg)),ones(size(origImg)), zeros(size(origImg))).*signal(:,:,sl);
red = cat(3, ones(size(origImg)),zeros(size(origImg)), zeros(size(origImg))).*interference(:,:,sl);
hold on
h1 = imshow(green); h2 = imshow(red);
set(h1, 'AlphaData', 0.15); set(h2, 'AlphaData', 0.15);
hold off
title('signal and interference regions');
set(gca, 'FontSize', 12)

% axial orientation
% figure(22);clf;
% sl = 65;
% origImg = flip(rot90(squeeze(img_orig_cc(sl,:,:)),1),2);
% imshow(origImg,[0 0.2]);
% green = cat(3, zeros(size(origImg)),ones(size(origImg)), zeros(size(origImg))).*flip(rot90(squeeze(signal(sl,:,:)),1),2);
% red = cat(3, ones(size(origImg)),zeros(size(origImg)), zeros(size(origImg))).*flip(rot90(squeeze(interference(sl,:,:)),1),2);
% hold on
% h1 = imshow(green); h2 = imshow(red);
% set(h1, 'AlphaData', 0.15); set(h2, 'AlphaData', 0.15);
% hold off
% title('signal and interference regions');
% set(gca, 'FontSize', 12)
% axis square

subplot 132
imshow(origImg,[0 0.2]);
title('all coils')
set(gca, 'FontSize', 12)

subplot 133
imshow(img_rovir_cc(:,:,sl),[0 0.2]);
title(sprintf('ROVir, %i coils',nv))
set(gca, 'FontSize', 12)
set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','FontName'),'FontName','Microsoft Yahei')
set(gcf,'Position',[964 53 954 913])
end

