function recon_dce_3Dsensemap( RP )
%RECONFUN Summary of this function goes here

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

% %% Eric, sorting the respiratory binning further for registration
% switch RP.respSortMethod
%     case 1 % use the RNAV
%         [extr2, nPhases] = retro_RNAV_RegistrationSort(RESP.Parameter.Labels, ...
%             RESP.Parameter.Cardiac.RNAV, RP.respPhases);
%     case 2  % use VitalEye
%         [extr2, nPhases] = retro_VitalEye_RegistrationSort(RESP.Parameter.Labels, ...
%             -double(RESP.Parameter.Labels.Index.na), RP.respPhases);
%     case 3  % use self-gating, takes time
%         logger.note(sprintf('Performing self-gating on centre k-space line'));  
%         %sg_signal_final = sg_rawSignal_aortaCINE(RESP);
%         sg_signal_final = sg_rawSignal(RESP);
%         [extr2, nPhases] = retro_SG_RegistrationSort(RESP.Parameter.Labels, ...
%             sg_signal_final, RP.respPhases);
%     case 4  % use PEAR belt signal, need the external scanphyslog file
%         pear_signal = readPEARSignal(RESP, RP);
%         [extr2, nPhases] = retro_PEAR_RegistrationSort(RESP.Parameter.Labels, ...
%             pear_signal, RP.respPhases);
% end
% 
% savename = fullfile(dir_out(RP), RP.name);
% mkdir(savename)
% imwrite(frame2im(getframe(figure(1))),[savename '/resp_bin_sorting.png']);
%% start recon

mrecon = RESP.Copy;
clear RESP
logger.note(sprintf('Read data'));
% read normal data only
mrecon.Parameter.Parameter2Read.typ = 1;
mrecon.Parameter.Parameter2Read.Update;
mrecon.ReadData;

% perform MRecon reconstruction steps
mrecon.RandomPhaseCorrection;
mrecon.PDACorrection;
mrecon.DcOffsetCorrection;
mrecon.MeasPhaseCorrection;

%% Perform sorting for DCE recon
logger.note(sprintf('Sorting into DCE frames'));
[extr2, nFrames] = sortDCEFrames(mrecon, RP.DCETimeResolution);

%%
% update parameters before sorting
mrecon.Parameter.Labels.Index.extr2 = extr2;
mrecon.Parameter.Parameter2Read.extr2 = (0:nFrames-1)';
mrecon.Parameter.Parameter2Read.Update;

mrecon.SortData;
mrecon.GridData;


%% 3D SENSEMAP
logger.note(sprintf('Calculate 3D sensemaps'));
tmp = mrecon.Data;

% undo checkerboard in mrecon.Data
tmp = bsxfun(@times,create_checkerboard([size(tmp,1),size(tmp,2),size(tmp,3)]),tmp);

% get sensitivity map
cmdsens = RP.bart_sensemap_cmd;
% average over all cardiac and resp bins
tmp3 = sum(sum(tmp(:,:,:,:,1,:,1,1,1,1,:),6),11)./sum(sum(tmp(:,:,:,:,1,:,1,1,1,1,:)~=0+eps,6),11);
tic
[L,sensemap] = bart_evalc(cmdsens, tmp3);
toc

shiftdata = fftshift(sensemap, 1);
sensemap = shiftdata;
clear tmp tmp3 shiftdata

% Crop the readout direction of the sensemaps
I = sensemap;
sz = size(I,1)/2;
y=(size(I,1)-sz)/2;
y=ceil(y+1:y+sz);
I=I(y,:,:,:);
sensemap = I;
clear I

% remove oversampling in kx and kz
mrecon.RemoveOversampling;

%%
% export sampling mask
logger.note(sprintf('Export sampling mask'));
try
    mask = mrecon.Data ~= 0;
    mask_dims = size(mask);
    pic = imstackmontage(squeeze(mask(mask_dims(1)/2,:,:,1,1,1,1,1,1,1,:,1)));
    m = mask(:,:,:,:,:,:,:,:,:,:,:,1);
    R = numel(m)/sum(m(:)) /4*pi;
    savename = fullfile(dir_out(RP), RP.name);
    imwrite(pic, sprintf('%s_mask_R%.1f.png',savename,R),'png');
    logger.note(sprintf('R by mask: %.2f', R));
    if checkRPflag(RP,'exp_mask') 
        save_mat(dir_out(RP), [RP.name '.mask'], 'mask', mask);
    end
    clear mask mask_dims pic m R
catch ME
    logger.note(sprintf('WARNING: Error occured in recon.m -- export sampling mask\n\n%s',ME.message'));
end

% place to reduce number of frames to speed up reconstruction
% mrecon.Data = mrecon.Data(:,:,:,:,:,:,:,:,:,:,1:20);

% pre-whitening data
if ~checkRPflag(RP,'reco_skip_PW')
    logger.note(sprintf('Pre-whitening data'));
    try
        MRn=MRecon(fullfile(RP.data_dir, RP.data_target));
        MRn.Parameter.Recon.ArrayCompression = mrecon.Parameter.Recon.ArrayCompression;
        MRn.Parameter.Recon.ACNrVirtualChannels = mrecon.Parameter.Recon.ACNrVirtualChannels;
        MRn.Parameter.Parameter2Read.typ=5;
        MRn.ReadData;
        eta=MRn.Data;

        Ncoils=size(mrecon.Data,4);
        Nsamples=numel(mrecon.Data)/Ncoils;

        psi = (1/(Nsamples-1))*(eta' * eta);
        L = chol(psi,'lower');
        L_inv = (inv(L));
        L_inv=diag(diag(L_inv)); %using only diagonal values

        mrecon.Data=permute(mrecon.Data, [1:3 5:length(size(mrecon.Data)) 4]);
        sizeMRDATA = size(mrecon.Data);
        mrecon.Data=reshape(mrecon.Data,[Nsamples,Ncoils]);
        mrecon.Data=mrecon.Data.';
        mrecon.Data = conj(L_inv) * mrecon.Data;
        mrecon.Data=mrecon.Data.';
        mrecon.Data=reshape(mrecon.Data,sizeMRDATA);
        mrecon.Data=ipermute(mrecon.Data, [1:3 5:length(size(mrecon.Data)) 4]);
        clearvars sizeMRDATA L_inv L psi MRn eta
    catch ME
            logger.note(sprintf('WARNING: Error occured in recon.m -- pre-whitening data\n\n%s',ME.message'));
    end
end

mrecon.RingingFilter;
% Eric, 20210104 check for k-space zerofilling, if 'No' saves memory
if isfield(RP,'Recon_kSpaceZeroFill')
    if strncmp(RP.Recon_kSpaceZeroFill, 'Yes',5)
        mrecon.ZeroFill;
    end
end

% coil compression / combination
if ( checkRPflag(RP,'bart_skip_fmac_SENSEUnfold') && checkRPflag(RP,'bart_cc') )
    try
        logger.note(sprintf('BART coil compression'));
        tmp = mrecon.Data;
        if isfield(RP,'bart_cc_cmd')
            cmdcc = RP.bart_cc_cmd;
        else
            cmdcc = 'cc -p 8 -E';
        end
        [T,tmp] = evalc('bart(cmdcc,tmp)');
        logger.note(T);
        mrecon.Data = tmp;
        clear tmp
    catch ME
        logger.note(sprintf('WARNING: Error occured in recon.m -- coil compression / combination!\n\n%s',ME.message'));
    end
elseif ( ~checkRPflag(RP,'bart_skip_fmac_SENSEUnfold') && checkRPflag(RP,'bart_cc') )
    logger.note(sprintf('INFO: BART coil compression cannot be performed if SENSEUnfold is enabled [ bart_skip_fmac_SENSEUnfold disabled ].'));
end


% load / estimate sensitivity maps
try
    if ( ~checkRPflag(RP,'bart_skip_fmac_SENSEUnfold') || ~checkRPflag(RP,'bart_sensemap') )
        % load MRsense object
        logger.note(sprintf('Load MRsense object'));
        S = MRsense(fullfile(RP.data_dir, RP.data_senseref), fullfile(RP.data_dir, RP.data_target));
        if checkRPflag(RP,'bart_pics')
            S.OutputSizeReformated = [size(mrecon.Data,1), ...
                                      size(mrecon.Data,2),...
                                      size(mrecon.Data,3)];
            S.OutputSizeSensitivity = S.OutputSizeReformated;
        end
        S.Mask = 1;
        S.Smooth = 1;
        S.Extrapolate = 1;
        S.Perform;
        mrecon.Parameter.Recon.Sensitivities = S;
        mrecon.Parameter.Recon.SENSERegStrength = 0;
        sensemap_all = mrecon.Parameter.Recon.Sensitivities.Sensitivity;
        if strcmp(mrecon.Parameter.Recon.ArrayCompression,'Yes')
            Sensitivity = S.Sensitivity;
            s_size = size(Sensitivity);
            Sensitivity = reshape(permute(Sensitivity,[4 1 2 3]),[s_size(4) prod(s_size(1:3))]);
            ACMatrix = mrecon.Parameter.Recon.ACMatrix;
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

%% perform CS reconstruction 
mrecon.K2IM;
try
    if ( checkRPflag(RP,'bart_sensemap') || checkRPflag(RP,'bart_pics') )      
        % save local tmp data
        logger.note(sprintf('Save k-space locally'));
        tmpsavename = ['/scratch/emschrauben' tempname];
        mkdir(tmpsavename)
        n_FE = size(mrecon.Data,1);
        for i_FE = 1:n_FE
            tmp = mrecon.Data(i_FE,:,:,:,:,:,:,:,:,:,:,:);
            tmpsavename_full = fullfile(tmpsavename,sprintf('slice_%03d',i_FE));
            writecfl(tmpsavename_full,tmp);
        end
        if ~checkRPflag(RP,'bart_skip_fmac_SENSEUnfold')
            tmp_out = mrecon.Data;
        else
            tmp_out = mrecon.Data(:,:,:,1,:,:,:,:,:,:,:,:);
        end
        if ~exist('sensemap_all','var'); sensemap_all = mrecon.Data(:,:,:,:,1,1,1,1,1,1,1,1); end
        mrecon.Data = [];
                
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
            tmpsavename_full = fullfile(tmpsavename,sprintf('slice_%03d',i_FE));
            tmp = readcfl(tmpsavename_full);
            
            % undo checkerboard in mrecon.Data
            tmp = bsxfun(@times,create_checkerboard([1,size(tmp,2),size(tmp,3)]),tmp);  
        
            % BART: PICS
            try
                if checkRPflag(RP,'bart_pics')
                    if isfield(RP,'bart_pics_cmd')
                        cmdpics = RP.bart_pics_cmd;
                    else
                        cmdpics = 'pics -R W:7:0:0.01 -R T:1024:0:0.01 -i 50 -S -d 5';
                    end
                    tmp = permute(tmp, dims_change_mrecon2bart);
                    % [BART MRI DIMS: READ_DIM,	PHS1_DIM,	PHS2_DIM,	COIL_DIM,	MAPS_DIM,	TE_DIM,	COEFF_DIM,	COEFF2_DIM,	ITER_DIM,	CSHIFT_DIM,	TIME_DIM,	TIME2_DIM,	LEVEL_DIM,	SLICE_DIM,	AVG_DIM
                    currSenseMap = sensemap(i_FE,:,:,:);
                    [L,tmp] = bart_evalc(cmdpics,tmp,currSenseMap);
                    if i_FE ==1 || i_FE == floor(n_FE/2); logger.note(L); end
                    if ~checkRPflag(RP,'bart_skip_fmac_SENSEUnfold')
                        tmp = bart('fft -u 7',bart('fmac -s 16',tmp,sensemap)); % first fmac than fft
                    end
                    tmp = ipermute(tmp, dims_change_mrecon2bart);
%                     tmp = bsxfun(@times,create_checkerboard([1,size(tmp,2),size(tmp,3)]),tmp); % redo checkerboard like in mrecon.Data
                    tmp_out(i_FE,:,:,:,:,:,:,:,:,:,:,:) = tmp;
                end
            catch ME
                logger.note(sprintf('WARNING: Error occured in recon.m -- BART: PICS!\n\n%s',ME.message'));
            end
        end % end of parfor single slice loop
        delete(gcp('nocreate'))

        % delete local tmp data
        for i_FE = 1:n_FE
            if (exist(sprintf(fullfile(tmpsavename,'/slice_%03d.hdr'),i_FE),'file'))
                delete(sprintf(fullfile(tmpsavename,'/slice_%03d.hdr'),i_FE));
            end
            if (exist(sprintf(fullfile(tmpsavename,'/slice_%03d.cfl'),i_FE),'file'))
                delete(sprintf(fullfile(tmpsavename,'/slice_%03d.cfl'),i_FE));
            end
        end
        rmdir(tmpsavename);
        mrecon.Data = tmp_out;
        sensemap_all = sensemap;
        clear tmp_out tmp sensemap T
    end
catch ME
    logger.note(sprintf('WARNING: Error occured in recon.m -- perform CS reconstruction\n\n%s',ME.message'));
end

if ( ~checkRPflag(RP,'bart_skip_fmac_SENSEUnfold') || ~checkRPflag(RP,'bart_pics') )
    logger.note(sprintf('SENSEUnfold correction'));
    mrecon.EPIPhaseCorrection;
    mrecon.K2IP;
    mrecon.GridderNormalization;
    mrecon.SENSEUnfold;
    
else
    try
        % CLEAR correction for CS reco
        if ~checkRPflag(RP,'reco_skip_CLEAR')
            logger.note(sprintf('CLEAR correction for CS reco'));
            mr_dims = size(mrecon.Data);
            data = squeeze(mrecon.Data);
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
            mrecon.Data = reshape(data_corr,mr_dims);
        end
        mrecon.Parameter.ReconFlags.isimspace  = [1,1,1];
        mrecon.Parameter.ReconFlags.isdepicorr = 1;
        mrecon.Parameter.ReconFlags.isunfolded = 1;
    catch ME
        logger.note(sprintf('WARNING: Error occured in recon.m -- CLEAR correction for CS reco!\n\n%s',ME.message'));
    end
end
%%
if ~checkRPflag(RP,'bart_pics'); mrecon.PartialFourier; end

mrecon.ConcomitantFieldCorrection;
mrecon.CombineCoils;
mrecon.Average;
mrecon.GeometryCorrection;

% manually remove outer slices here?

% mrecon.RemoveOversampling;
% mrecon.ZeroFill;
mrecon.FlowPhaseCorrection;
mrecon.RotateImage;

% Rescale image data
mrecon.Data = data_rescale(mrecon.Data);

%% Clip noise
try
    if isfield(RP,'reco_NoiseClip') && strcmp(RP.reco_NoiseClip, 'OnlyPCMRA')
        if ~isfield(RP,'reco_NoiseClipValue'); RP.reco_NoiseClipValue = 2^8; end
        NCmrecon = mrecon.Copy;
        NCmrecon.Data = data_mask(mrecon.Data,RP.reco_NoiseClipValue);
    end
catch ME
    logger.note(sprintf('WARNING: Error occured in recon.m -- Clip noise!\n\n%s',ME.message'));
end


%% export data
logger.note(sprintf('Export data'));
savename = fullfile(dir_out(RP), RP.name);

% save MRecon flags
try
    param_list.Scan = properties(mrecon.Parameter.Scan);
    param_list.Recon = properties(mrecon.Parameter.Recon);
    fieldlist = fieldnames(param_list);
    fileID = fopen(fullfile(dir_out(RP), [RP.name '.MReconFlags.dat']), 'w');
    for ii = 1 : size(fieldlist,1)
        strtable = eval(sprintf('param_list.%s',fieldlist{ii,1}));
        for jj = 1 : eval(sprintf('size(param_list.%s,1);',fieldlist{ii,1}))
            if eval(sprintf('isa(mrecon.Parameter.%s.%s,''double'')',fieldlist{ii},strtable{jj,1}))
                fprintf(fileID,'%-35s\t%s\n',[fieldlist{ii} '.' strtable{jj,1}], eval(sprintf('num2str(mrecon.Parameter.%s.%s)',fieldlist{ii},strtable{jj,1} )) );
            end
            if eval(sprintf('isa(mrecon.Parameter.%s.%s,''char'')',fieldlist{ii},strtable{jj,1}))
                fprintf(fileID,'%-35s\t%s\n',[fieldlist{ii} '.' strtable{jj,1}], eval(sprintf('mrecon.Parameter.%s.%s',fieldlist{ii},strtable{jj,1} )) );
            end
        end
    end
    fclose(fileID);
    clear param_list fieldlist ii jj fileID strtable ans
catch ME
    logger.note(sprintf('WARNING: Error occured in recon.m -- save MRecon flags -- data could not be saved!\n%s',ME.message'));
end

% save .png image
try
    plotphase = 1;
    magnitude = abs(mrecon.Data(:,:,:,1,1,1,1,1,1,1,plotphase,1,1));
    imwrite(imstackmontage(magnitude) ./ prctile(magnitude(:),95), [savename,'_Mag_M',num2str(plotphase),'.png'],'png'); %scaled to 95 percentile
catch ME
    logger.note(sprintf('WARNING: Error occured in recon.m -- save .png image -- data could not be saved!\n%s',ME.message'));
end

% save sensitivity maps
if checkRPflag(RP,'exp_sense')
    try
        save_mat( fullfile(dir_out(RP)), [RP.name '.Sensitivity'], 'Sensitivity',sensemap_all);
    catch ME
        logger.note(sprintf('WARNING: Error occured in recon.m -- save sensitivity maps -- data could not be saved!\n%s',ME.message'));
    end
end

% export loop over respiratory dimension
data_storage = mrecon.Data;
for i_resp = 1:size(mrecon.Data,11)
    mrecon.Data = data_storage(:,:,:,:,:,:,:,:,:,:,i_resp,:);
    % save .par/.rec
    if checkRPflag(RP,'exp_parrec')
        try
            mrecon.WritePar([savename, '_resp', num2str(i_resp), '.par']);
            mrecon.WriteRec([savename, '_resp', num2str(i_resp), '.rec']);
            if exist('pcmra','var') && i_resp==1
                pcmra.WritePar([savename '_PCMRA.par']);
                pcmra.WriteRec([savename '_PCMRA.rec']);
            end
        catch ME
            logger.note(sprintf('WARNING: Error occured in recon.m -- save .par/.rec -- data could not be saved!\n%s',ME.message'));
        end
    end
    
    % save .dcm
    if checkRPflag(RP,'exp_dcm')
        try
            % manual dicom export to fix missing dicom tags
            
            % create Parfiles and OutputDirectories
            parcell = mrecon.WritePar2String;
            reccell = mrecon.GetRecImages;
            npar = length(parcell);
            mpscell = {'M', 'P', 'S'};
            for ipar = 1:npar
                currentParfileCell{ipar} = mrecon.Parameter.parread_from_string(parcell{ipar});
                OutputDirectoryCell{ipar} = [dir_out(RP),'/DICOM_',mpscell{ipar}]; % ,'_resp',num2str(i_resp)]; removed by Carmen
                rawFileCell{ipar} = mrecon.Parameter.Par40;
            end
            if exist('pcmra','var') && i_resp==1
                npar = npar+1;
                OutputDirectoryCell{npar} = fullfile(dir_out(RP),'PCMRA');
                parcellpcmra = pcmra.WritePar2String;
                parcell(npar) = parcellpcmra;
                currentParfileCell{npar} = pcmra.Parameter.parread_from_string(parcell{npar});
                reccell(npar) = pcmra.GetRecImages;
                rawFileCell{npar} = pcmra.Parameter.Par40;
            end
            % perform DICOMExporterUser in parallel
            if npar > 1; parpool(npar); end
            parfor ipar = 1:npar;
                D = DICOMExporter;
                mkdir(OutputDirectoryCell{ipar});
                D.Export(OutputDirectoryCell{ipar}, reccell{ipar}, currentParfileCell{ipar}, rawFileCell{ipar});
            end
            delete(gcp('nocreate'))
        catch ME
            logger.note(sprintf('WARNING: Error occured in recon.m -- save .dcm -- data could not be saved!\n%s',ME.message));
        end
    end
end
mrecon.Data = data_storage;
clear data_store

% save .mat
if checkRPflag(RP,'exp_mat')
    try
        mrecon.Data=squeeze(mrecon.Data);
        mrecon.Data(isnan(mrecon.Data)) = 0;
        save_mat( fullfile(dir_out(RP)), [RP.name '.Data'], 'Data',mrecon.Data);
    catch ME
        logger.note(sprintf('WARNING: Error occured in recon.m -- save .mat -- data could not be saved!\n%s',ME.message'));
    end
end

% save .nii
if checkRPflag(RP,'exp_nii')
    %lmg debug: this part seems not to work
    try
        addpath('/home/emschrauben/lood_storage/divi/Pub/Software/matlab/toolbox_common/nifti/20140122/')
        magnitude = squeeze(abs(mrecon.Data)); % lmg debug: is this rescaling needed?
        magnitude = magnitude ./ max(magnitude(:));
        magnitude = make_nii(magnitude, mrecon.Parameter.Scan.RecVoxelSize);
        save_nii(magnitude,[savename '-mag.nii']);
    catch ME
        logger.note(sprintf('WARNING: Error occured in recon.m -- save .nii -- data could not be saved!\n%s',ME.message'));
    end
end

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