%% Load resp stuff, redo registration
clc; clear all; close all

addpath(genpath(fullfile(cd,'./lib/')));
[fName, pName] = uigetfile(['/home/emschrauben/lood_storage/divi/Projects/asap/*.mat']);
load(fullfile(pName,fName))
%%
nPhases = size(magnitude,4);
% 3D registration of only the pancreas
recon = squeeze(mean(magnitude,4));
recon = recon/max(recon(:));
temp = squeeze(recon(:,:,round(size(recon,3)/2)));
figure(1); clf; imagesc(temp); colormap gray, axis equal off
h = drawrectangle(gca);
pause;
maskSz = round(h.Position);
close(figure(1))
%     save([RP.data_dir 'RespMask.mat'],'maskSz')

% take only the central 50% of slices
magnitude2 = magnitude(maskSz(2):(maskSz(2)+maskSz(4)), maskSz(1):(maskSz(1)+maskSz(3)),...
    round(0.25*size(magnitude,3)):round(0.75*size(magnitude,3)),:);

clear transforms registered
for phase = 2:nPhases
    fixed = magnitude2(:,:,:,1);
    moving = magnitude2(:,:,:,phase);
    [T, RegisteredImage] = registerImages3D(moving,fixed);
    transforms(:,phase-1) = squeeze(T.T(4,1:3))';    % only store the translational components
    registered(:,:,:,phase-1) = RegisteredImage;
end

% to view the results
a = repmat(fixed,[1 1 1 nPhases]);
b = cat(4,fixed,registered);
View4D(cat(2,magnitude2,a,b),1,'axisnames',...
    {'','','moving, end-expiration, registered'}, 'FigureName','Registration Result',...
    'FramePanelTitle','Resp Frames')

% save([savename '/transforms.mat'],'transforms');
% save([savename '/respRecon.mat'],'magnitude','magnitude2','registered','extr2','maskSz');

transforms = round(transforms,2);
tra_table = table((1:size(transforms,2))',squeeze(transforms(1,:))',squeeze(transforms(2,:))',squeeze(transforms(3,:))',...
    'VariableNames',{'phase','AP','RL','FH'});
tra_table
% writetable(tra_table,[savename '/translations.xlsx']);