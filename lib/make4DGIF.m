clc; clear all; close all

if isunix
    addpath('/home/emschrauben/lood_storage/divi/Pub/Software/matlab/toolbox_common/nifti/20140122/')
    addpath(genpath('/home/emschrauben/scratch/code/matlab/toolboxes/'))
    [fName, pName] = uigetfile(['/home/emschrauben/lood_storage/divi/Projects/asap/*.nii']);
else
    addpath('L:\basic\divi\Pub\Software\matlab\toolbox_common\nifti\20140122\')
    addpath('L:\basic\divi\Projects\dcerecon\Spirals_for_DCE\matlab code\rscaspr-master\Simulations\autofocus\View4D')
    [fName, pName] = uigetfile(['L:\basic\divi\Projects\asap\*.nii']);
end


tmp = load_nii(fullfile(pName,fName));

mag = tmp.img; mag = mag/max(abs(mag(:)));

View4D(abs(squeeze(flip(mag,3))),1,'PixelDimensions',tmp.hdr.dime.pixdim(2:4),...
    'IntersectCoordinates',[58,123,15])