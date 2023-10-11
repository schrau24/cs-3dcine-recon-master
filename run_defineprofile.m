% run defineprofile

% copy default profile (defineprofile_default.m) to defineprofile.m and open it


% initialization
addpath(fullfile(cd,'./lib/'));
clear;
close all;
clc;

% copy and open
copyfile(fullfile(cd,'lib','defineprofile_default.m'),fullfile(cd,'lib','defineprofile.m'));
open('defineprofile.m');