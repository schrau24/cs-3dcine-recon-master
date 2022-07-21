function [extr2, transforms] = recon_RespRegister( RP )
%RECONFUN Summary of this function goes here
%

%% presets
% check data file type
if ~strcmp( RP.data_target(end-3:end) , '.raw') && ~strcmp( RP.data_target(end-3:end) , '.lab')
    RP.data_target = strcat(RP.data_target,'.raw');
end
if ~strcmp( RP.data_senseref(end-3:end) , '.raw') && ~strcmp( RP.data_senseref(end-3:end) , '.lab')
    RP.data_senseref = strcat(RP.data_senseref,'.raw');
end
if ~strcmp( RP.data_coilsurvey(end-3:end) , '.raw') && ~strcmp( RP.data_coilsurvey(end-3:end) , '.lab')
    RP.data_coilsurvey = strcat(RP.data_coilsurvey,'.raw');
end

% initialize MRecon object
init_toolbox;
global logger
RESP = MRecon(fullfile(RP.data_dir, RP.data_target)); % '.raw' or '.lab'

% set flags to disable/enable certain MRrecon operations
RPlist = fields(RP);
for ii=2:size(RPlist,1)
    strcell = strsplit(RPlist{ii,1},'_');
    if strcmp(strcell{1,1},'Recon') || strcmp(strcell{1,1},'Cardiac') || strcmp(strcell{1,1},'Encoding')
        eval(sprintf('RESP.Parameter.%s.%s = RP.%s;',strcell{1,1},strcell{1,2},RPlist{ii,1}))
    end
end
clear RPlist strcell

%% extract respiratory self-gating signal
if checkRPflag(RP,'respCorrection')
    logger.note(sprintf('Performing self-gating on centre k-space line'));
    [RespSignal, sg_signal_final] = respiratory_motion(RESP, RP.DCETimeResolution > 30);
    [extr2, nPhases] = retro_SG_Sort(RESP.Parameter.Labels,sg_signal_final,RP.respPhases);
    savename = fullfile(dir_out(RP), RP.name);
    mkdir(savename)
    imwrite(frame2im(getframe(figure(1))),[savename '/resp_signal.png']);
end

%% Read in data
logger.note(sprintf('Read data'));

% read normal data only
RESP.Parameter.Parameter2Read.typ = 1;
RESP.Parameter.Parameter2Read.Update;
RESP.ReadData;

% perform MRecon reconstruction steps
RESP.RandomPhaseCorrection;
RESP.RemoveOversampling;
RESP.PDACorrection;
RESP.DcOffsetCorrection;
RESP.MeasPhaseCorrection;

%%

% update parameters before sorting
RESP.Parameter.Labels.Index.extr2 = extr2;
RESP.Parameter.Parameter2Read.extr2 = (0:nPhases-1)';
RESP.Parameter.Parameter2Read.Update;

RESP.SortData;
RESP.GridData;

% export sampling mask
logger.note(sprintf('Export sampling mask'));
try
    mask = RESP.Data ~= 0;
    mask_dims = size(mask);
    pic = imstackmontage(squeeze(mask(mask_dims(1)/2,:,:,1,1,1,1,1,1,1,:,1)));
    m = mask(:,:,:,:,:,:,:,:,:,:,:,1);
    R = numel(m)/sum(m(:)) /4*pi;
    savename = fullfile(dir_out(RP), RP.name);
    imwrite(pic, sprintf('%s_respmask_R%.1f.png',savename,R),'png');
    logger.note(sprintf('R by mask: %.2f', R));
    if checkRPflag(RP,'exp_mask') 
        save_mat(dir_out(RP), [RP.name '.mask'], 'mask', mask);
    end
    clear mask mask_dims pic m R
catch ME
    logger.note(sprintf('WARNING: Error occured in recon.m -- export sampling mask\n\n%s',ME.message'));
end

% pre-whitening data
if ~checkRPflag(RP,'reco_skip_PW')
    logger.note(sprintf('Pre-whitening data'));
    try
        MRn=MRecon(fullfile(RP.data_dir, RP.data_target));
        MRn.Parameter.Recon.ArrayCompression = RESP.Parameter.Recon.ArrayCompression;
        MRn.Parameter.Recon.ACNrVirtualChannels = RESP.Parameter.Recon.ACNrVirtualChannels;
        MRn.Parameter.Parameter2Read.typ=5;
        MRn.ReadData;
        eta=MRn.Data;
        
        Ncoils=size(RESP.Data,4);
        Nsamples=numel(RESP.Data)/Ncoils;
        
        psi = (1/(Nsamples-1))*(eta' * eta);
        L = chol(psi,'lower');
        L_inv = (inv(L));
        L_inv=diag(diag(L_inv)); %using only diagonal values
        
        RESP.Data=permute(RESP.Data, [1:3 5:length(size(RESP.Data)) 4]);
        sizeMRDATA = size(RESP.Data);
        RESP.Data=reshape(RESP.Data,[Nsamples,Ncoils]);
        RESP.Data=RESP.Data.';
        RESP.Data = conj(L_inv) * RESP.Data;
        RESP.Data=RESP.Data.';
        RESP.Data=reshape(RESP.Data,sizeMRDATA);
        RESP.Data=ipermute(RESP.Data, [1:3 5:length(size(RESP.Data)) 4]);
        clearvars sizeMRDATA L_inv L psi MRn eta
    catch ME
        logger.note(sprintf('WARNING: Error occured in recon.m -- pre-whitening data\n\n%s',ME.message'));
    end
end

RESP.RingingFilter;
% Eric, 20210104 check for k-space zerofilling, if 'No' saves memory
if isfield(RP,'Recon_kSpaceZeroFill')
    if strncmp(RP.Recon_kSpaceZeroFill, 'Yes',5)
        RESP.ZeroFill;
    end
end

% load / estimate sensitivity maps from senserefscan
try
    if ~checkRPflag(RP,'bart_sensemap')
        % load MRsense object
        logger.note(sprintf('Load MRsense object'));
        S = MRsense(fullfile(RP.data_dir, RP.data_senseref), fullfile(RP.data_dir, RP.data_target));
        if checkRPflag(RP,'bart_pics')
            S.OutputSizeReformated = [size(RESP.Data,1), ...
                                      size(RESP.Data,2),...
                                      size(RESP.Data,3)];
            S.OutputSizeSensitivity = S.OutputSizeReformated;
        end
        S.Mask = 1;
        S.Smooth = 1;
        S.Extrapolate = 1;
        S.Perform;
        RESP.Parameter.Recon.Sensitivities = S;
        RESP.Parameter.Recon.SENSERegStrength = 0;
        sensemap_all = RESP.Parameter.Recon.Sensitivities.Sensitivity;
        if strcmp(RESP.Parameter.Recon.ArrayCompression,'Yes')
            Sensitivity = S.Sensitivity;
            s_size = size(Sensitivity);
            Sensitivity = reshape(permute(Sensitivity,[4 1 2 3]),[s_size(4) prod(s_size(1:3))]);
            ACMatrix = RESP.Parameter.Recon.ACMatrix;
            ac_size = size(ACMatrix);
            SensitivityAC = ACMatrix(:,1:s_size(4)) * Sensitivity;
            SensitivityAC = permute(reshape(SensitivityAC,[ac_size(1) s_size(1:3)]),[2 3 4 1]);
            sensemap_all = SensitivityAC;
            S.Sensitivity = SensitivityAC;
            clear Sensitivity s_size ACMatrix ac_size SensitivityAC
        end        
        if ~checkRPflag(RP,'bart_sensemap')
            sensemap_all = sensemap_all./max(sensemap_all(:)); % LMG DEBUG (normalize for better image quality)
        end
        clear S
    end
catch ME
    logger.note(sprintf('WARNING: Error occured in recon.m -- load MRsense object!\n\n%s',ME.message'));  
end   

RESP.K2IM;
%% perform CS reconstruction
try
    if ( checkRPflag(RP,'bart_sensemap') || checkRPflag(RP,'bart_pics') )
        % save local tmp data
        logger.note(sprintf('Save k-space locally'));
        tmpsavename = tempname;
        mkdir(tmpsavename)
        n_FE = size(RESP.Data,1);
        for i_FE = 1:n_FE
            tmp = RESP.Data(i_FE,:,:,:,:,:,:,:,:,:,:,:);
            save_mat(tmpsavename,sprintf('slice_%03d',i_FE),'tmp',tmp);
        end
        if ~checkRPflag(RP,'bart_skip_fmac_SENSEUnfold')
            tmp_out = RESP.Data;
        else
            tmp_out = RESP.Data(:,:,:,1,:,:,:,:,:,:,:,:);
        end
        if ~exist('sensemap_all','var'); sensemap_all = RESP.Data(:,:,:,:,1,1,1,1,1,1,1,1); end
        RESP.Data = [];
        
        % create parpool
        delete(gcp('nocreate'))
        if ~isfield(RP,'bart_parpool'); RP.bart_parpool = 8; end
        if RP.bart_parpool > 0
            logger.note(sprintf('Create parallel pool (%d workers)',RP.bart_parpool));
            parpool(RP.bart_parpool);
        else
            ps = parallel.Settings;
            ps.Pool.AutoCreate =  false;
        end
        
        logger.note(sprintf('CS reconstruction'));
        parfor i_FE = 1:n_FE % start of parfor single slice loop
            display(sprintf('CS recon: slice %03d',i_FE));
            % load local tmp slice
            T = load(sprintf(fullfile(tmpsavename,'/slice_%03d.mat'),i_FE));
            tmp = T.tmp;
            
            % undo checkerboard in RESP.Data
            tmp = bsxfun(@times,create_checkerboard([1,size(tmp,2),size(tmp,3)]),tmp);
            
            try % get sensitivity map
                if checkRPflag(RP,'bart_sensemap')
                    % BART estimate sensitivity
                    if isfield(RP,'bart_sensemap_cmd')
                        cmdsens = RP.bart_sensemap_cmd;
                    else
                        cmdsens = 'ecalib -m1';
                    end
                    [L,sensemap] = bart_evalc(cmdsens, sum(sum(sum(tmp(:,:,:,:,1,:,1,1,1,:),6),10),11)./sum(sum(sum(tmp(:,:,:,:,1,:,1,1,1,:)~=0+eps,6),10),11) );
                    sensemap_all(i_FE,:,:,:) = sensemap;
                    if i_FE ==1 || i_FE == floor(n_FE/2); logger.note(L); end
                else
                    % load MRSense map
                    sensemap = sensemap_all(i_FE,:,:,:);
                end
            catch ME
                logger.note(sprintf('WARNING: Error occured in recon.m -- BART estimate sensitivity!\n\n%s',ME.message'));
            end
            
            % BART: PICS
            try
                if checkRPflag(RP,'bart_pics')
                    
                    % manually input BART pics command
                    % respiratory temporal dimension
                    cmdpics = 'pics -R T:6:0:0.01 -R T:2048:0:0.01 -i 50 -S -d 5';
                    tmp = permute(tmp, dims_change_mrecon2bart);
                    % [BART MRI DIMS: READ_DIM,	PHS1_DIM,	PHS2_DIM,	COIL_DIM,	MAPS_DIM,	TE_DIM,	COEFF_DIM,	COEFF2_DIM,	ITER_DIM,	CSHIFT_DIM,	TIME_DIM,	TIME2_DIM,	LEVEL_DIM,	SLICE_DIM,	AVG_DIM
                    
                    [L,tmp] = bart_evalc(cmdpics,tmp,sensemap);
                    if i_FE ==1 || i_FE == floor(n_FE/2); logger.note(L); end
                    if ~checkRPflag(RP,'bart_skip_fmac_SENSEUnfold')
                        tmp = bart('fft -u 7',bart('fmac -s 16',tmp,sensemap)); % first fmac than fft
                    end
                    %                     tmp = permute(single(tmp), [1:5 11 6:9 10]); % permute tmp back for RESP.Data
                    %                     tmp = permute(single(tmp), [1:4 7 11 5 6 8 9 10]); % permute tmp back for RESP.Data (working rnav implementation in dyn dimension)
                    tmp = ipermute(tmp, dims_change_mrecon2bart);
                    tmp = bsxfun(@times,create_checkerboard([1,size(tmp,2),size(tmp,3)]),tmp); % redo checkerboard like in RESP.Data
                    tmp_out(i_FE,:,:,:,:,:,:,:,:,:,:,:) = tmp;
                end
            catch ME
                logger.note(sprintf('WARNING: Error occured in recon.m -- BART: PICS!\n\n%s',ME.message'));
            end
        end % end of parfor single slice loop
        delete(gcp('nocreate'))
        
        % delete local tmp data
        for i_FE = 1:n_FE
            delete(sprintf(fullfile(tmpsavename,'/slice_%03d.mat'),i_FE));
        end
        rmdir(tmpsavename);
        RESP.Data = tmp_out;
        clear tmp_out tmp sensemap T
    end
catch ME
    logger.note(sprintf('WARNING: Error occured in recon.m -- perform CS reconstruction\n\n%s',ME.message'));
end
try
    % CLEAR correction for CS reco
    if ~checkRPflag(RP,'reco_skip_CLEAR')
        logger.note(sprintf('CLEAR correction for CS reco'));
        mr_dims = size(RESP.Data);
        data = squeeze(RESP.Data);
        data_dims = size(data);
        l2norm=sqrt(sum(abs(sensemap_all).^2,4));
        l2norm_repmat = repmat(l2norm,[1 1 1 data_dims(4:end)]);
        l2norm_repmat_mask = l2norm_repmat<1.0;
        l2norm_repmat(l2norm_repmat_mask)=1;
        data_corr = data./(l2norm_repmat.^2);
        if isfield(RP,'reco_NoiseClip') && (strcmp(RP.reco_NoiseClip, 'All'));
            if ~isfield(RP,'reco_NoiseClipValue');
                RP.reco_NoiseClipValue = 2^8;
            end
            data_M_mask = abs(data_mask(data_rescale(data),RP.reco_NoiseClipValue)) > 0;
            data_corr = data_corr .* data_M_mask;
        end % this if-statement was added by Carmen
        RESP.Data = reshape(data_corr,mr_dims);
    end
    RESP.Parameter.ReconFlags.isimspace  = [1,1,1];
    RESP.Parameter.ReconFlags.isdepicorr = 1;
    RESP.Parameter.ReconFlags.isunfolded = 1;
catch ME
    logger.note(sprintf('WARNING: Error occured in recon.m -- CLEAR correction for CS reco!\n\n%s',ME.message'));
end

%% Calculate PC MRA
rr = RESP.Copy;
rr.ConcomitantFieldCorrection;
rr.DivideFlowSegments;
% rr.GeometryCorrection;
rr.FlowPhaseCorrection;

magnitude  =   abs(squeeze(mean(rr.Data(:,:,:,1,1,1,1,1,1,:,:,1),10)));
magnitude = rescale(magnitude);

% 3D registration of only the pancreas
recon = squeeze(mean(magnitude,4));
recon = recon/max(recon(:));

%% Register respiratory phases together

if isfile([RP.data_dir 'RespMask.mat'])
    load([RP.data_dir 'RespMask.mat']); % see if mask already exists
else
    temp = squeeze(recon(:,:,round(size(recon,3)/2)));
    figure(1); clf; imagesc(temp); colormap gray, axis equal off
    h = drawrectangle(gca);
    pause;
    maskSz = round(h.Position);
    close(figure(1))
    save([RP.data_dir 'RespMask.mat'],'maskSz')
end

% take only the central 50% of slices
magnitude2 = magnitude(maskSz(2):(maskSz(2)+maskSz(4)), maskSz(1):(maskSz(1)+maskSz(3)),...
    round(0.1*size(magnitude,3)):round(0.9*size(magnitude,3)),:);

clear transforms registered
for phase = 2:size(magnitude,4)
    fixed = magnitude2(:,:,:,1);
    moving = magnitude2(:,:,:,phase);
    [T, RegisteredImage] = registerImages3D(moving,fixed);
    transforms(:,phase-1) = squeeze(T.T(4,1:3))';    % only store the translational components
    registered(:,:,:,phase-1) = RegisteredImage;
end

% % to view the results
% a = repmat(fixed,[1 1 1 nPhases]);
% b = cat(4,fixed,registered);
% View4D(cat(2,magnitude2,a,b),1,'axisnames',...
%     {'','','moving, end-expiration, registered'}, 'FigureName','Registration Result',...
%     'FramePanelTitle','Resp Frames')

save([savename '/transforms.mat'],'transforms');
save([savename '/respRecon.mat'],'magnitude','magnitude2','registered','extr2','maskSz');

transforms = round(transforms,2);
tra_table = table((1:size(transforms,2))',squeeze(transforms(1,:))',squeeze(transforms(2,:))',squeeze(transforms(3,:))',...
    'VariableNames',{'phase','AP','RL','FH'});
tra_table
writetable(tra_table,[savename '/translations.xlsx']);
end

function B = isnat( x )
%ISNAT Summary of this function goes here
%   Detailed explanation goes here

B = (x == floor(x)) && (x > 0);
end

function out = checkRPflag( RP, field )
% check if profile parameter exists and is true
% 'field' needs to be a string

out = isfield(RP,field);
if out
    eval(sprintf('out = (RP.%s == 1);',field));
end
end

function ch=create_checkerboard(s)
%s: size of checkerboard
% starts with -1 on top left corner
if length(s)==2
    ch=(((-1).^[1:s(1)]).*1i).'*(((-1).^[1:s(2)]).*1i);
elseif length(s)==3
    ch=(((-1).^[1:s(1)]).*1i).'*(((-1).^[1:s(2)]).*1i);
    ch=repmat(ch,[1 1 s(3)]);
    ch1d=(((-1).^[1:s(3)]).*1i).*(ones(1,s(3)).*-1i);
    ch1d=permute(ch1d,[1 3 2]);
    ch=bsxfun(@times,ch,ch1d);
elseif length(s)==1
    ch=(((-1).^[1:s(1)]).*1i).*(ones(1,s(1)).*-1i);
else
    error('unsupported size')
end

end

function dims_order = dims_change_mrecon2bart()
% change mrecon dimension to bart dimensions
%
% MRECON MRI DIMS: x  y  z  coils  dynamics cardiac_phases  echoes  locations  mixes  extr1  extr2  averages
% BART MRI DIMS: READ_DIM,	PHS1_DIM,	PHS2_DIM,	COIL_DIM,	MAPS_DIM,	TE_DIM,	COEFF_DIM,	COEFF2_DIM,	ITER_DIM,	CSHIFT_DIM,	TIME_DIM,	TIME2_DIM,	LEVEL_DIM,	SLICE_DIM,	AVG_DIM

mrecon2bart = [ 1, 1;...% READ_DIM      (x)
    2, 2; ...           % PHS1_DIM      (y)
    3, 3; ...           % PHS2_DIM      (z)
    4, 4; ...           % COIL_DIM      (coils)
    5, 8; ...           % MAPS_DIM      (locations)
    6, 7; ...           % TE_DIM        (echoes)
    7, 5; ...           % COEFF_DIM     (dynamics)
    8, 9; ...           % COEFF2_DIM    (mixes)
    9, 10; ...          % ITER_DIM      (extr1)
    10, 13; ...         % CSHIFT_DIM    (NONE)
    11, 6; ...          % TIME_DIM      (cardiac_phases)
    12, 11; ...         % TIME2_DIM     (extr2)
    13, 14; ...         % LEVEL_DIM     (NONE)
    14, 15; ...         % SLICE_DIM     (NONE)
    15, 12]';           % AVG_DIM       (averages)

dims_order = mrecon2bart(2,:);
end

function data = data_rescale(data)
% out = data_rescale(in)
%
% scale data to 2^12 = 4096

datamax = ceil(max(abs(data(:))));
data = data ./ datamax .* 2^12;
end

function data = data_mask(data,NoiseClipValue)
% data = data_mask(data)
%
% mask the data based on a noise level

% take mean over entire data set (x,y,z, : )
dims = size(data);
th_data = mean(abs(reshape(data,[dims(1:3),prod(dims(4:end))])),4);

% make mask by threshold / NoiseClipValue
th_mask = th_data > NoiseClipValue;

% apply mask to data
th_mask_all = repmat(th_mask,[1,1,1,dims(4:end)]);
data = data .* th_mask_all;
end