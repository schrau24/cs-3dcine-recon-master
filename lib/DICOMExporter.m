classdef DICOMExporter < handle
    
    properties
        
    end
    
    properties (Hidden, SetAccess = protected )
        mDictionary = [];
        mInfoStatic = struct;
        mInfo = struct;
        mGoalcPars = [];
        mRelease = 'REL3';
        mOutDir = [];
        mOutDirRelativePath = [];
        mSequences = [];
        mPatientID = [];         
    end
    
    methods
        function D = DICOMExporter()
            dicomdict('set', 'MRdict.txt');
            D.define_sequences();
        end
        
        function SetDictionary( D, dic )
            try
                D.mDictionary = dic;
                dicomdict('set', dic);
            catch
                error( 'Cannot set the DICOM dictionary' );
            end
        end
        
        function Export(D, OutDir, Recfile, Parfile, GoalcPars)
            if( nargin < 4 )
                error( 'Not enough input arguments');
            end
            
            if( ~iscell(Parfile) )
                Parfile = {Parfile};
            end
            if( ~iscell(Recfile) )
                Recfile = {Recfile};
            end
            
            if( length(Recfile) ~= length(Parfile))
                error( 'The same number of parameter- and datasets must be given as input');
            end
            
            if( nargin == 5 && ~strcmpi(class(GoalcPars), 'ParameterReader'))
                error( 'The 5th argument must be a ParfileReader class');
            end
            
            if( ~isdir( OutDir ))
                error( [OutDir, ' is not a directory'] );
            end
                        
            OutDir = strrep(OutDir, '/', filesep);
            OutDir = strrep(OutDir, '\', filesep);
            
            D.mOutDir = OutDir;
            
            % LMG edit: not needed
%             if( ~isdir( [OutDir, filesep, 'DICOM'] ))
%                 try
%                     mkdir(OutDir, 'DICOM');
%                 catch
%                     error( ['cannot create the directory OutDir ', OutDir, filesep, 'DICOM'] );
%                 end
%             end           
            
            
            if( nargin == 5 )
                D.mGoalcPars = GoalcPars;
                if (D.mGoalcPars.IsParameter('RFR_EXAM_OBJECT_OID'))
                    if( D.mGoalcPars.IsParameter('VAL01_ACQ_echoes'))
                        D.mRelease = 'REL53';
                    else
                        D.mRelease = 'REL5';
                    end
                end
            end
            
            D.get_filenames(length(Parfile));
            
            % create a patient id
            time_format = 'HHMMSS.FFF';
            cur_time = datestr(now,time_format);
            D.mPatientID =  [Parfile{1}.PatientName(1:2), '_',  cur_time];
            
            for i = 1:length(Parfile)
                
                par = Parfile{i};
                if( ischar(par) )
                    % the parfile is a filename
                    ParfileStruct = MRparameter.parread(par);
                elseif (isstruct(par ))
                    % the parfile is already a struct
                    ParfileStruct = par;
                end
                rec = Recfile{i};
                if( ischar(rec) )
                    Data = MRecon.readrec(rec);
                elseif (isnumeric(rec ))
                    % the parfile is already a struct
                    Data = rec;
                end
                
                
                % Delete the info structs
                D.mInfo = struct;
                D.mInfoStatic = struct;
                
                fprintf( 'DICOMExporter: Getting parameters...');
                if( ~isempty(D.mGoalcPars) ) 
                    if strcmp(D.mRelease, 'REL53')
                        D.dicom_get_static_values53(ParfileStruct);
                    else
                        D.dicom_get_static_values(ParfileStruct);
                    end
                else
                    D.dicom_get_static_values_from_parfile(ParfileStruct);
                end
                D.dicom_write_images(Data, ParfileStruct, i);
            end
        end
    end
    
    methods (Hidden, Access = protected )
        function dicom_write_images(D, Data, Parfile, file_number)
            if( ~isempty(D.mGoalcPars) )
                try
                    if( strcmpi( D.mRelease(1:4), 'REL5'))
                        series_uid = D.mGoalcPars.GetValue('RFR_SERIES_DICOM_SERIES_INSTANCE_UID');
                    else
                        series_uid = D.mGoalcPars.GetValue('RFR_SERIES_series_instance_uid');
                    end
                    series_uid = D.get_uid( series_uid );
                catch
                    series_uid = dicomuid;
                end
            else
                series_uid = dicomuid;
            end
            
            cur_series_uid = series_uid;
            cur_series_uid(end) = num2str(file_number);
            
            nr_images = size(Data(:,:,:), 3);
            if( length(Parfile.ImageInformation.SliceNumber) ~= nr_images)
                error( 'The number of images in the data and the parfile do not match');
            end
            
            counter = Counter( 'DICOMExporter: Writing image %d/%d');
            nr_images = size(Data(:,:,:), 3);
            
            % LMG DEBUG: test filter
            
            
            
            % LMG DEBUG start
            % sort Data into Philips order
            id_philips = [];
            [list_ImageTypeMr,~]=sort(Parfile.ImageInformation.ImageTypeMr);
            unique_ImageTypeMr=unique(list_ImageTypeMr);
            for ii = 1:length(unique_ImageTypeMr)
                [~,idx_SliceNumber]=sort(Parfile.ImageInformation.SliceNumber(list_ImageTypeMr==unique_ImageTypeMr(ii)));
                id_philips = [id_philips;idx_SliceNumber+length(id_philips)];
            end
            Parfile.ImageInformation.IndexInRECFile(id_philips) = Parfile.ImageInformation.IndexInRECFile;
            Data = Data(:,:,id_philips);
            % LMG DEBUG end
            
            for image_nr = 1:nr_images
                D.mInfo = D.mInfoStatic;
                
                cur_type = Parfile.ImageInformation.ImageTypeMr(image_nr);
                
                if( ~isempty(D.mGoalcPars) )
                    D.dicom_get_dynamic_values(Parfile, image_nr, cur_type, file_number)
                else
                    D.dicom_get_dynamic_values_from_parfile(Parfile, image_nr)
                end
                
%                 cur_dir = [D.mOutDir, filesep, 'DICOM']; %LMG edit
                % LMG edit: not needed
                cur_dir = D.mOutDir;
                
                if( file_number <= length(D.mOutDirRelativePath) && ~isempty(D.mOutDirRelativePath{file_number} ) )
                    cur_dir = [cur_dir, filesep, D.mOutDirRelativePath{file_number} ];
                end
                
                % create the directories
                if( ~isdir( cur_dir ) )
                    parent_dir = [D.mOutDir, filesep, 'DICOM'];
                    mkdir(parent_dir, D.mOutDirRelativePath{file_number});
                end
                
                % LMG DEBUG: filename as InstanceNumber / IndexInRECFile
                image_nr_export = Parfile.ImageInformation.IndexInRECFile(image_nr)+1;
                if (image_nr_export < 10)
                    temp_file = ['IM_000', num2str(image_nr_export)];
                elseif (image_nr_export < 100)
                    temp_file = ['IM_00', num2str(image_nr_export)];
                elseif (image_nr_export < 1000)
                    temp_file = ['IM_0', num2str(image_nr_export)];
                elseif (image_nr_export < 10000)
                    temp_file = ['IM_', num2str(image_nr_export)];
                end
                
                cur_filename = [cur_dir, filesep, temp_file ];
                
                D.mInfo = DICOMExporter.sort_tags(D.mInfo);
                
                % Update Counter
                counter.Update( {image_nr, nr_images} );    
                % use the save method when writing release 3 DICOM's
                % (Matlab will verify the tags)
                dicomwrite(Data(:,:,image_nr_export)', cur_filename, D.mInfo, 'CreateMode', 'copy', 'WritePrivate', true)
            end           
            DICOMExporter.create_dicomdir( D.mOutDir, D.mOutDir );
            fprintf('\n');
        end
        function dicom_get_static_values(D, Parfile)
            D.dicom_get_static_values_from_parfile(Parfile);
            D.dicom_get_database_values();
            
            try
                % Get the rest of the static values
                scan_sequence = 'RFR_SERIES_PIIM_SERIES_SCAN_SEQUENCE';
                sequence_name = 'RFR_SERIES_PIIM_MR_SERIES_SCANNING_TECHNIQUE_DESC';
                field_strength = 'RFR_SERIES_PIIM_MR_SERIES_MAGNETIC_FIELD_STRENGTH';
                scan_duration = 'RFR_SERIES_PIIM_MR_SERIES_SCAN_DURATION';
                epi_factor = 'RFR_SERIES_PIIM_MR_SERIES_EPI_FACTOR';
                phase_contrast_nr_dirs = 'RFR_SERIES_PIIM_MR_SERIES_NR_OF_PHASE_CONTRAST_DIRCTNS';
                fat_saturation = 'RFR_SERIES_PIIM_MR_FAT_SATURATION';
                
                if (strcmpi(D.mRelease, 'REL3'))
                    TR = [D.mGoalcPars.GetValue('RFR_SERIES_repetition_times0'); 0];
                    series_uid = D.get_uid(D.mGoalcPars.GetValue('RFR_SERIES_series_instance_uid'));
                    birth_date = format_date(D.mGoalcPars.GetValue('RFR_PATIENT_patient_birth_date'));
                    scan_sequence = 'RFR_SERIES_series_scan_sequence';
                    sequence_name = 'RFR_SERIES_pulse_sequence_name';
                    field_strength = 'RFR_SERIES_magnetic_field';
                    scan_duration = 'RFR_SERIES_scan_duration';
                    epi_factor = 'RFR_SERIES_epi_factor';
                    phase_contrast_nr_dirs = 'RFR_SERIES_no_pc_directions';
                    fat_saturation = 'RFR_SERIES_fat_saturation';
                    acquisition_type_num = D.mGoalcPars.GetValue('RFR_SERIES_mr_acquisition_type');
                    
                    switch (str2double(acquisition_type_num))
                        case 1
                            acquisition_type = '2D';
                        case 2
                            acquisition_type = '3D';
                        case 3
                            acquisition_type = 'MS';
                        case 4
                            acquisition_type = 'M2D';
                        case 5
                            acquisition_type = 'SV';
                        case 6
                            acquisition_type = '1D';
                    end
                else                    
                    TR = [str2double(D.mGoalcPars.GetValue('RFR_SERIES_PIIM_MR_SERIES_REPETITION_TIME')); 0];
                    if isnan(TR(1))
                        try
                            TR = [D.mGoalcPars.GetValue('RFR_SERIES_PIIM_MR_SERIES_REPETITION_TIME'); 0];
                        catch end
                    end
                    series_uid = D.get_uid(D.mGoalcPars.GetValue('RFR_SERIES_DICOM_SERIES_INSTANCE_UID'));
                    birth_date = D.mGoalcPars.GetValue('RFR_PATIENT_DICOM_PATIENT_BIRTH_DATE');
                    acquisition_type = D.mGoalcPars.GetValue('RFR_SERIES_PIIM_MR_SERIES_ACQUISITION_TYPE_PRIVATE');
                end
                
                vCurPar = 'PatientAge';  D.set_par(vCurPar, DICOMExporter.get_patients_age(birth_date),'static');
                
                vCurPar = 'SeriesInstanceUID';  D.set_par(vCurPar, series_uid,'static');
                vCurPar = {'0x2005', '0x1030'};  D.set_par(vCurPar, TR,'static');
                
                vCurPar = 'CodingSchemeDesignator';  D.set_par(vCurPar, 'DCM','static');
                vCurPar = 'SamplesPerPixel';  D.set_par(vCurPar, 1,'static');
                vCurPar = 'VolumetricProperties';  D.set_par(vCurPar, 'VOLUME','static');
                vCurPar = 'VolumeBasedCalculationTechnique';  D.set_par(vCurPar, 'NONE','static');
                vCurPar = 'TimeDomainFiltering';  D.set_par(vCurPar, '','static');
                vCurPar = 'MetaboliteMapDescription';  D.set_par(vCurPar, 'WATER','static');
                vCurPar = 'KSpaceFiltering';  D.set_par(vCurPar, 'HAMMING','static');
                vCurPar = 'LUTLabel';  D.set_par(vCurPar, 'Philips','static');
                vCurPar = {'0x2005', '0x1035'};  D.set_par(vCurPar, 'PIXEL','static'); 																																																																																																% still empty (no CRecon parameter found)
                vCurPar = 'ReferringPhysicianName'; D.set_par(vCurPar, '','static');
                vCurPar = 'AdditionalPatientHistory'; D.set_par(vCurPar, '','static');
                vCurPar = 'SequenceVariant'; D.set_par(vCurPar, '','static');
                vCurPar = 'ScanOptions'; D.set_par(vCurPar, '','static');
                vCurPar = 'MultiCoilElementName'; D.set_par(vCurPar, '','static');
                vCurPar = 'DecoupledNucleus'; D.set_par(vCurPar, '','static');
                
                echo_pulse_sequence = 'GRADIENT';
                steady_state_pulse_seq = 'NONE';
                
                if ( ~isempty(strfind(D.mGoalcPars.GetValue(sequence_name), 'MIX')))
                    echo_pulse_sequence = 'MIXED';
                else
                    if (strcmpi(D.mGoalcPars.GetValue(scan_sequence),'GR') || ...
                            strcmpi(D.mGoalcPars.GetValue(scan_sequence),'1') )
                        echo_pulse_sequence = 'GRADIENT';
                        if ( ~isempty(strfind(D.mGoalcPars.GetValue(sequence_name),'T1')))
                            steady_state_pulse_seq = 'LONGITUDINAL';
                        elseif ( ~isempty(strfind(D.mGoalcPars.GetValue(sequence_name),'T2')))
                            steady_state_pulse_seq = 'TIME_REVERSE';
                        elseif ( ~isempty(strfind(D.mGoalcPars.GetValue(sequence_name),'B-')))
                            steady_state_pulse_seq = 'FREE_PRECESSION';
                        end
                    elseif (strcmpi(D.mGoalcPars.GetValue(scan_sequence),'SE') || ...
                            strcmpi(D.mGoalcPars.GetValue(scan_sequence),'0') )
                        echo_pulse_sequence = 'SPIN';
                    end
                end
                vCurPar = 'EchoPulseSequence';  D.set_par(vCurPar, echo_pulse_sequence,'static');
                vCurPar = 'SteadyStatePulseSequence';  D.set_par(vCurPar, steady_state_pulse_seq,'static');
                
                if( strcmpi(acquisition_type, '3D'))
                    vCurPar = 'MRAcquisitionType';  D.set_par(vCurPar, '3D','static');
                else
                    vCurPar = 'MRAcquisitionType';  D.set_par(vCurPar, '2D','static');
                end
                
                vCurPar = 'MagneticFieldStrength';  D.set_par(vCurPar, D.mGoalcPars.GetValue(field_strength),'static');
                vCurPar = 'AcquisitionDuration';  D.set_par(vCurPar, D.mGoalcPars.GetValue(scan_duration),'static');
                
                % SpectrallySelectedExcitation
                if (strcmpi(D.mGoalcPars.GetValue(fat_saturation),'N') || ...
                        strcmpi(D.mGoalcPars.GetValue(fat_saturation),'0') )
                    SpectrallySelectedExcitation = 'NONE';
                else
                    SpectrallySelectedExcitation = 'FAT';
                end
                vCurPar = 'SpectrallySelectedExcitation';  D.set_par(vCurPar, SpectrallySelectedExcitation,'static');
                
                if( str2double(D.mGoalcPars.GetValue(phase_contrast_nr_dirs)) > 1 )
                    vCurPar = 'PhaseContrast';  D.set_par(vCurPar, 'YES','static');
                else
                    vCurPar = 'PhaseContrast';  D.set_par(vCurPar, 'NO','static');
                end
                if( str2double(D.mGoalcPars.GetValue(epi_factor)) > 1 )
                    vCurPar = 'EchoPlanarPulseSequence';  D.set_par(vCurPar, 'YES','static');
                else
                    vCurPar = 'EchoPlanarPulseSequence';  D.set_par(vCurPar, 'NO','static');
                end
                
                if (D.mGoalcPars.IsParameter('EX_ACQ_imaging_sequence'))
                    
                    switch (D.mGoalcPars.GetValue('UGN1_FFE_contrast_enh', [], 1))
                        case 0 %MGUFFE_CE_NO:
                            steady_state_pulse_seq = 'NONE';
                        case 1 % MGUFFE_CE_T1:
                            steady_state_pulse_seq = 'LONGITUDINAL';
                        case 2 % MGUFFE_CE_T2:
                            steady_state_pulse_seq = 'TIME_REVERSE';
                        case 3 % MGUFFE_CE_HYBRID:
                            steady_state_pulse_seq = 'FREE_PRECESSION';
                        otherwise
                            steady_state_pulse_seq = 'NONE';
                    end
                    
                    % Tagging
                    if (D.mGoalcPars.GetValue('UGN1_TAG_ctagging') ~= 0)
                        
                        switch (D.mGoalcPars.GetValue('UGN1_TAG_dimension', [], 1))
                            
                            case 0 % MGU_TAG_DIM_LINES:
                                tagging = 'LINE';
                            case 1 %MGU_TAG_DIM_GRID:
                                tagging = 'GRID';
                            case 2 %MGU_TAG_DIM_2XLINES:
                                tagging = 'GRID';
                            otherwise
                                tagging = 'NONE';
                        end
                    else
                        tagging = 'NONE';
                    end
                    
                    % MTC
                    switch (D.mGoalcPars.GetValue('UGN1_MTC_enable', [], 1))
                        case 0 % MPU_MTC_MODE_NO:
                            mtc = 'NONE';
                        case 1 %MPU_MTC_MODE_ON_RES:
                            mtc = 'ON_RESONANCE';
                        case 2 %MPU_MTC_MODE_OFF_RES:
                            mtc = 'OFF_RESONANCE';
                        case 3 %MPU_MTC_MODE_OFF_RES_MULTI_PULSE:
                            mtc = 'OFF_RESONANCE';
                        otherwise
                            mtc = 'NONE';
                    end
                    
                    % SpectrallySelectedExcitation
                    if (D.mGoalcPars.GetValue('UGN1_SPIR_mode') == 3) % MPUSPIR_MODE_PROSET)
                        SpectrallySelectedExcitation = 'NONE';
                    else
                        if (D.mGoalcPars.GetValue('UGN1_SPIR_suppr_tissue') == 0)
                            SpectrallySelectedExcitation = 'WATER';
                        else
                            SpectrallySelectedExcitation = 'FAT';
                        end
                    end
                    
                    % SpectrallySelectedSuppression
                    SpectrallySelectedSuppression = 'NONE';
                    water_suppress = false;
                    fat_suppress = false;
                    if (D.mGoalcPars.GetValue('UGN1_MTC_sos') == 1)
                        water_suppress = true;
                        fat_suppress = true;
                    end
                    
                    if (D.mGoalcPars.GetValue('UGN1_SPIR_enable') == 1)
                        if (D.mGoalcPars.GetValue('UGN1_SPIR_enable') == 1) %MPUSPIR_SUPPR_WATER
                            water_suppress = true;
                        else
                            fat_suppress = true;
                        end
                    end
                    
                    if (D.mGoalcPars.GetValue('UGN1_WATSUP_method') ~= 0)% MGUWATSUP_METHOD_NO
                        water_suppress = true;
                    end
                    
                    switch (D.mGoalcPars.GetValue('UGN1_WATSUP_in_base', [], 1))
                        case 2 % MPUWATSUP_INBASE_WATER:
                            water_suppress = true;
                        case 1 % MPUWATSUP_INBASE_FAT:
                            fat_suppress = true;
                        case 3 % MPUWATSUP_INBASE_WATER_FAT:
                            water_suppress = true;
                            fat_suppress = true;
                    end
                    
                    if (water_suppress && fat_suppress)
                        SpectrallySelectedSuppression = 'FAT_AND_WATER';
                    end
                    
                    if (D.mGoalcPars.GetValue('UGN1_WATSUP_method') ~= 0) %MGUWATSUP_METHOD_NO
                        water_suppress = true;
                    end
                    
                    switch (D.mGoalcPars.GetValue('UGN1_WATSUP_in_base', [], 1))
                        case 2 % MPUWATSUP_INBASE_WATER:
                            water_suppress = true;
                        case 1 % MPUWATSUP_INBASE_FAT:
                            fat_suppress = true;
                        case 3 %MPUWATSUP_INBASE_WATER_FAT:
                            water_suppress = true;
                            fat_suppress = true;
                    end
                    
                    if (water_suppress && ~fat_suppress)
                        SpectrallySelectedSuppression = 'WATER';
                    end
                    if (~water_suppress && fat_suppress)
                        SpectrallySelectedSuppression = 'FAT';
                    end
                    if (water_suppress && fat_suppress)
                        SpectrallySelectedSuppression = 'FAT_AND_WATER';
                    end
                    
                    % GeometryOfKSpaceTraversal
                    switch (D.mGoalcPars.GetValue('RC_k_space_traj_type', [], 1))
                        case 0
                            GeometryOfKSpaceTraversal = 'RECTILINEAR';
                        case 1
                            GeometryOfKSpaceTraversal = 'RADIAL';
                        case 2
                            GeometryOfKSpaceTraversal = 'SPIRAL';
                        case 3
                            GeometryOfKSpaceTraversal = 'PROPELLER';
                        otherwise
                            GeometryOfKSpaceTraversal = 'RECTILINEAR';
                    end
                    
                    % RectilinearPhaseEncodeReordering
                    if (D.mGoalcPars.GetValue('UGN1_ACQ_imaging_sequence') == 3) % MGUACQ_SEQ_FFE
                        outer_order = D.mGoalcPars.GetValue('UGN1_FFE_outer_profile_orders', 1, 1);
                    else
                        outer_order = D.mGoalcPars.GetValue('UGN1_RFE_outer_profile_orders', 1, 1);
                    end
                    
                    switch (outer_order)
                        case 2 % MGOBJORD_PROF_ORD_LINEAR:
                            RectilinearPhaseEncodeReordering = 'LINEAR';
                        case 3 % MGOBJORD_PROF_ORD_REV_LINEAR:
                            RectilinearPhaseEncodeReordering = 'REVERSE_LINEAR';
                        case 4 % MGOBJORD_PROF_ORD_SWEEP:
                            RectilinearPhaseEncodeReordering = 'SEGMENTED';
                        case 8 % MGOBJORD_PROF_ORD_SWEEP_2:
                            RectilinearPhaseEncodeReordering = 'SEGMENTED';
                        case 9 % MGOBJORD_PROF_ORD_SWEEP_4:
                            RectilinearPhaseEncodeReordering = 'SEGMENTED';
                        case 5 % MGOBJORD_PROF_ORD_OPT_CENTRIC:
                            RectilinearPhaseEncodeReordering = 'CENTRIC';
                        case 7 % MGOBJORD_PROF_ORD_OPT_MIN_CENTRIC:
                            RectilinearPhaseEncodeReordering = 'CENTRIC';
                        case 1 % MGOBJORD_PROF_ORD_LOW_HIGH:
                            RectilinearPhaseEncodeReordering = 'CENTRIC';
                        case 6 % MGOBJORD_PROF_ORD_OPT_REV_CENTRIC:
                            RectilinearPhaseEncodeReordering = 'REVERSE_CENTRIC';
                        case 0 % MGOBJORD_PROF_ORD_HIGH_LOW:
                            RectilinearPhaseEncodeReordering = 'REVERSE_CENTRIC';
                        case 10 % MGOBJORD_PROF_ORD_UNDEF:
                            RectilinearPhaseEncodeReordering = 'UNKNOWN';
                        case 12 % MGOBJORD_PROF_ORD_RADIAL:
                            RectilinearPhaseEncodeReordering = 'UNKNOWN';
                        case 13 % MGOBJORD_PROF_ORD_SPIRAL:
                            RectilinearPhaseEncodeReordering = 'UNKNOWN';
                        case 14 % MGOBJORD_PROF_ORD_VISTA:
                            RectilinearPhaseEncodeReordering = 'UNKNOWN';
                        otherwise
                            RectilinearPhaseEncodeReordering = 'UNKNOWN';
                    end
                    
                    % ParallelAcquisitionTechnique
                    ParallelAcquisitionTechnique = '';
                    if (D.mGoalcPars.GetValue('UGN1_GEO_sense_enable', [], 1) == 1)
                        ParallelAcquisitionTechnique = 'SENSE';
                    end
                    if (D.mGoalcPars.GetValue('UGN1_ACQ_kt_factor') > 1)
                        ParallelAcquisitionTechnique = 'KT';
                    end
                    if (D.mGoalcPars.IsParameter('UGN1_ACQ_MB') && (D.mGoalcPars.GetValue('UGN1_ACQ_MB', [], 1) == 1))
                        ParallelAcquisitionTechnique = 'MULTIBAND';
                    end
                    if (D.mGoalcPars.IsParameter('UGN1_ACQ_caipi') && (D.mGoalcPars.GetValue('UGN1_ACQ_caipi', [], 1) == 1))
                        ParallelAcquisitionTechnique = 'MULTIBAND';
                    end
                    
                    % CardiacSignalSource
                    CardiacSignalSource = '';
                    if (D.mGoalcPars.GetValue('UGN1_CARD_synchronisation', [], 1) ~= 0)
                        
                        switch (D.mGoalcPars.GetValue('EX_CARD_device', [], 1))
                            case 0 % MGG_CARD_DEV_ECG
                                CardiacSignalSource = 'ECG';
                            case 2 % MGG_CARD_DEV_INT
                                CardiacSignalSource = 'ECG';
                            case 3 % MGG_CARD_DEV_EXT
                                CardiacSignalSource = 'ECG';
                            case 1 % MGG_CARD_DEV_PPU
                                CardiacSignalSource = 'PP';
                        end
                    end
                    
                    % CoverageOfKSpace
                    if (D.mGoalcPars.GetValue('UGN1_ACQ_scan_mode', [], 1) == 1)
                        is_shutter = D.mGoalcPars.GetValue('UGN0_DEF_elliptical_k_space_shutter', [], 1) == 1;
                        if (is_shutter)
                            if (D.mGoalcPars.GetValue('UGN1_ACQ_scan_type', [], 1) == 1) % MGUACQ_SCT_SPECTRO
                                CoverageOfKSpace = 'ELLIPSOIDAL';
                            else
                                CoverageOfKSpace = 'ELLIPTICAL';
                            end
                            
                        else
                            CoverageOfKSpace = 'FULL';
                        end
                    else
                        CoverageOfKSpace = '';
                    end
                    
                    if (D.mGoalcPars.IsParameter('UGN1_DIFF_nr_directions'))
                        NrGradOris = D.mGoalcPars.GetValue('UGN1_DIFF_nr_directions');
                    elseif (D.mGoalcPars.IsParameter('UGN5_DIFF_measured_nr_oris'))
                        NrGradOris = D.mGoalcPars.GetValue('UGN5_DIFF_measured_nr_oris');
                    else
                        NrGradOris = 1;
                    end
                    
                    % PatientPosition
                    if( D.mGoalcPars.GetValue('EX_GEO_patient_position', [], 1) == 0)
                        patient_position = 'HF';
                    else
                        patient_position = 'FF';
                    end
                    if( D.mGoalcPars.GetValue('EX_GEO_patient_orientation', [], 1) == 0)
                        patient_position = [patient_position, 'S'];
                    elseif( D.mGoalcPars.GetValue('EX_GEO_patient_orientation', [], 1) == 1)
                        patient_position = [patient_position, 'P'];
                    elseif( D.mGoalcPars.GetValue('EX_GEO_patient_orientation', [], 1) == 2)
                        patient_position = [patient_position, 'L'];
                    else
                        patient_position = [patient_position, 'R'];
                    end
                    
                    vCurPar = 'PatientPosition'; D.set_par(vCurPar, patient_position ,'static');
                    
                    if( D.mGoalcPars.GetValue('UGN1_ACQ_scan_mode', [], 1) == 1 )
                        vCurPar = 'MRAcquisitionType';  D.set_par(vCurPar, '3D','static');
                    else
                        vCurPar = 'MRAcquisitionType';  D.set_par(vCurPar, '2D','static');
                    end
                    
                    vCurPar = 'MagneticFieldStrength';  D.set_par(vCurPar, D.mGoalcPars.GetValue('HW_main_magnetic_field_mT') / 1000,'static');
                    vCurPar = 'PercentPhaseFieldOfView';  D.set_par(vCurPar, D.mGoalcPars.GetValue('EX_GEO_rect_fov_perc'),'static');
                    vCurPar = 'AcquisitionDuration';  D.set_par(vCurPar, D.mGoalcPars.GetValue('AC_total_scan_time'),'static');
                    vCurPar = 'AcquisitionContrast';  D.set_par(vCurPar, D.mGoalcPars.GetValue('EX_FFE_ceffe'),'static');
                    
                    if( D.mGoalcPars.GetValue('UGN1_ACQ_fast_imaging_mode', [], 1) == 1 )
                        vCurPar = 'MultipleSpinEcho';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'MultipleSpinEcho';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    if( D.mGoalcPars.IsParameter('UGN1_ACQ_MB') && D.mGoalcPars.GetValue('UGN1_ACQ_MB', [], 1) == 1 )
                        vCurPar = 'MultiPlanarExcitation';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'MultiPlanarExcitation';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    if( D.mGoalcPars.GetValue('EX_PC_angio_mode', [], 1) == 2 )
                        vCurPar = 'PhaseContrast';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'PhaseContrast';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    if( D.mGoalcPars.GetValue('EX_PC_angio_mode', [], 1) == 1 )
                        vCurPar = 'TimeOfFlightContrast';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'TimeOfFlightContrast';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    if( D.mGoalcPars.GetValue('MP_RFE_rf_spoil_mode', [], 1) == 1 && D.mGoalcPars.GetValue('UGN1_FFE_contrast_enh', [], 1) == 1 )
                        vCurPar = 'Spoiling';  D.set_par(vCurPar, 'RF','static');
                    else
                        vCurPar = 'Spoiling';  D.set_par(vCurPar, 'NONE','static');
                    end
                    
                    vCurPar = 'SteadyStatePulseSequence';  D.set_par(vCurPar, steady_state_pulse_seq,'static');
                    vCurPar = 'Tagging';  D.set_par(vCurPar, tagging,'static');
                    
                    if( D.mGoalcPars.GetValue('UGN1_ACQ_fast_imaging_mode', [], 1) == 3 || D.mGoalcPars.GetValue('UGN1_ACQ_fast_imaging_mode', [], 1) == 5 || D.mGoalcPars.GetValue('UGN1_ACQ_fast_imaging_mode', [], 1) == 4 )
                        vCurPar = 'EchoPlanarPulseSequence';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'EchoPlanarPulseSequence';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    vCurPar = 'MagnetizationTransfer';  D.set_par(vCurPar, mtc,'static');
                    
                    if( D.mGoalcPars.GetValue('UGN1_T2PREP_enable', [], 1) == 1 )
                        vCurPar = 'T2Preparation';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'T2Preparation';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    if( D.mGoalcPars.GetValue('UGN1_TFEPP_type', [], 1) == 3 )
                        vCurPar = 'BloodSignalNulling';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'BloodSignalNulling';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    if( D.mGoalcPars.GetValue('UGN1_ACQ_tfe_excitations', [], 1) == 1 && D.mGoalcPars.GetValue('UGN1_TFEPP_type', [], 1) == 2 )
                        vCurPar = 'SaturationRecovery';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'SaturationRecovery';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    vCurPar = 'SpectrallySelectedExcitation';  D.set_par(vCurPar, SpectrallySelectedExcitation,'static');
                    vCurPar = 'SpectrallySelectedSuppression';  D.set_par(vCurPar, SpectrallySelectedSuppression,'static');
                    vCurPar = 'GeometryOfKSpaceTraversal';  D.set_par(vCurPar, GeometryOfKSpaceTraversal,'static');
                    
                    if( D.mGoalcPars.GetValue('RC_quad_receive_coils', [], 1) ~= 0 )
                        vCurPar = 'QuadratureReceiveCoil';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'QuadratureReceiveCoil';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    vCurPar = 'RectilinearPhaseEncodeReordering';  D.set_par(vCurPar, RectilinearPhaseEncodeReordering,'static');
                    
                    if( D.mGoalcPars.GetValue('RC_receive_coil_types', [], 1) == 3 )
                        vCurPar = 'ReceiveCoilType';  D.set_par(vCurPar, 'MULTICOIL','static');
                    else
                        if( D.mGoalcPars.GetValue('RC_receive_coil_types', [], 1) == 4 )
                            vCurPar = 'ReceiveCoilType';  D.set_par(vCurPar, 'UNKNOWN','static');
                        else
                            vCurPar = 'ReceiveCoilType';  D.set_par(vCurPar, strrep(D.mGoalcPars.GetValue('RC_receive_coil_types', 1), 'ARGRC_COIL_TYPE_', ''),'static');
                        end
                    end
                    
                    if( D.mGoalcPars.GetValue('UGN1_SPY_decouple', [], 1) ~= 0 )
                        vCurPar = 'Decoupling';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'Decoupling';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    vCurPar = 'ParallelAcquisitionTechnique';  D.set_par(vCurPar, ParallelAcquisitionTechnique,'static');
                    vCurPar = 'CardiacSignalSource';  D.set_par(vCurPar, CardiacSignalSource,'static');
                    vCurPar = 'CoverageOfKSpace';  D.set_par(vCurPar, CoverageOfKSpace,'static');
                    vCurPar = {'0x2001', '0x1009'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('VW_prepulse_delay'),'static');
                    vCurPar = {'0x2005', '0x1415'};  D.set_par(vCurPar, NrGradOris,'static');
                    
                    vCurPar = {'0x2005', '0x1054'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_ap_angulations', 1,'static'));
                    vCurPar = {'0x2005', '0x1055'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_fh_angulations', 1,'static'));
                    vCurPar = {'0x2005', '0x1056'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_rl_angulations', 1,'static'));
                    vCurPar = {'0x2005', '0x1057'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_ap_fovs', 1,'static'));
                    vCurPar = {'0x2005', '0x1058'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_fh_fovs', 1,'static'));
                    vCurPar = {'0x2005', '0x1059'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_rl_fovs', 1,'static'));
                    vCurPar = {'0x2005', '0x105a'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_ap_offcentres', 1,'static'));
                    vCurPar = {'0x2005', '0x105b'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_fh_offcentres', 1,'static'));
                    vCurPar = {'0x2005', '0x105c'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_rl_offcentres', 1,'static'));
                    vCurPar = {'0x2005', '0x105e'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_slab_orientations', 1,'static'));
                    vCurPar = {'0x2005', '0x105d'};  D.set_par(vCurPar, upper(D.mGoalcPars.GetValue('PS_volume_geom_types', 1,'static')));
                    vCurPar = {'0x2001', '0x1001'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('RC_chemical_shifts', 1,'static'));
                end
            catch
                warning('Could not get all DICOM tags. Please make sure you have the most recent version of the ReconFrame patch');
            end
        end
        function dicom_get_static_values53(D, Parfile)
            D.dicom_get_static_values_from_parfile(Parfile);
            D.dicom_get_database_values();
            
            try
                % Get the rest of the static values
                scan_sequence = 'RFR_SERIES_PIIM_SERIES_SCAN_SEQUENCE';
                sequence_name = 'RFR_SERIES_PIIM_MR_SERIES_SCANNING_TECHNIQUE_DESC';
                field_strength = 'RFR_SERIES_PIIM_MR_SERIES_MAGNETIC_FIELD_STRENGTH';
                scan_duration = 'RFR_SERIES_PIIM_MR_SERIES_SCAN_DURATION';
                epi_factor = 'RFR_SERIES_PIIM_MR_SERIES_EPI_FACTOR';
                phase_contrast_nr_dirs = 'RFR_SERIES_PIIM_MR_SERIES_NR_OF_PHASE_CONTRAST_DIRCTNS';
                fat_saturation = 'RFR_SERIES_PIIM_MR_FAT_SATURATION';
                
                if (strcmpi(D.mRelease, 'REL3'))
                    TR = [D.mGoalcPars.GetValue('RFR_SERIES_repetition_times0'); 0];
                    series_uid = D.get_uid(D.mGoalcPars.GetValue('RFR_SERIES_series_instance_uid'));
                    birth_date = format_date(D.mGoalcPars.GetValue('RFR_PATIENT_patient_birth_date'));
                    scan_sequence = 'RFR_SERIES_series_scan_sequence';
                    sequence_name = 'RFR_SERIES_pulse_sequence_name';
                    field_strength = 'RFR_SERIES_magnetic_field';
                    scan_duration = 'RFR_SERIES_scan_duration';
                    epi_factor = 'RFR_SERIES_epi_factor';
                    phase_contrast_nr_dirs = 'RFR_SERIES_no_pc_directions';
                    fat_saturation = 'RFR_SERIES_fat_saturation';
                    acquisition_type_num = D.mGoalcPars.GetValue('RFR_SERIES_mr_acquisition_type');
                    
                    switch (str2double(acquisition_type_num))
                        case 1
                            acquisition_type = '2D';
                        case 2
                            acquisition_type = '3D';
                        case 3
                            acquisition_type = 'MS';
                        case 4
                            acquisition_type = 'M2D';
                        case 5
                            acquisition_type = 'SV';
                        case 6
                            acquisition_type = '1D';
                    end
                else
                    TR = [str2double(D.mGoalcPars.GetValue('RFR_SERIES_PIIM_MR_SERIES_REPETITION_TIME')); 0];
                    series_uid = D.get_uid(D.mGoalcPars.GetValue('RFR_SERIES_DICOM_SERIES_INSTANCE_UID'));
                    birth_date = D.mGoalcPars.GetValue('RFR_PATIENT_DICOM_PATIENT_BIRTH_DATE');
                    acquisition_type = D.mGoalcPars.GetValue('RFR_SERIES_PIIM_MR_SERIES_ACQUISITION_TYPE_PRIVATE');
                end
                
                vCurPar = 'PatientAge';  D.set_par(vCurPar, DICOMExporter.get_patients_age(birth_date),'static');
                
                vCurPar = 'SeriesInstanceUID';  D.set_par(vCurPar, series_uid,'static');
                vCurPar = {'0x2005', '0x1030'};  D.set_par(vCurPar, TR,'static');
                
                vCurPar = 'CodingSchemeDesignator';  D.set_par(vCurPar, 'DCM','static');
                vCurPar = 'SamplesPerPixel';  D.set_par(vCurPar, 1,'static');
                vCurPar = 'VolumetricProperties';  D.set_par(vCurPar, 'VOLUME','static');
                vCurPar = 'VolumeBasedCalculationTechnique';  D.set_par(vCurPar, 'NONE','static');
                vCurPar = 'TimeDomainFiltering';  D.set_par(vCurPar, '','static');
                vCurPar = 'MetaboliteMapDescription';  D.set_par(vCurPar, 'WATER','static');
                vCurPar = 'KSpaceFiltering';  D.set_par(vCurPar, 'HAMMING','static');
                vCurPar = 'LUTLabel';  D.set_par(vCurPar, 'Philips','static');
                vCurPar = {'0x2005', '0x1035'};  D.set_par(vCurPar, 'PIXEL','static'); 																																																																																																% still empty (no CRecon parameter found)
                vCurPar = 'ReferringPhysicianName'; D.set_par(vCurPar, '','static');
                vCurPar = 'AdditionalPatientHistory'; D.set_par(vCurPar, '','static');
                vCurPar = 'SequenceVariant'; D.set_par(vCurPar, '','static');
                vCurPar = 'ScanOptions'; D.set_par(vCurPar, '','static');
                vCurPar = 'MultiCoilElementName'; D.set_par(vCurPar, '','static');
                vCurPar = 'DecoupledNucleus'; D.set_par(vCurPar, '','static');
                
                echo_pulse_sequence = 'GRADIENT';
                steady_state_pulse_seq = 'NONE';
                
                if ( ~isempty(strfind(D.mGoalcPars.GetValue(sequence_name), 'MIX')))
                    echo_pulse_sequence = 'MIXED';
                else
                    if (strcmpi(D.mGoalcPars.GetValue(scan_sequence),'GR') || ...
                            strcmpi(D.mGoalcPars.GetValue(scan_sequence),'1') )
                        echo_pulse_sequence = 'GRADIENT';
                        if ( ~isempty(strfind(D.mGoalcPars.GetValue(sequence_name),'T1')))
                            steady_state_pulse_seq = 'LONGITUDINAL';
                        elseif ( ~isempty(strfind(D.mGoalcPars.GetValue(sequence_name),'T2')))
                            steady_state_pulse_seq = 'TIME_REVERSE';
                        elseif ( ~isempty(strfind(D.mGoalcPars.GetValue(sequence_name),'B-')))
                            steady_state_pulse_seq = 'FREE_PRECESSION';
                        end
                    elseif (strcmpi(D.mGoalcPars.GetValue(scan_sequence),'SE') || ...
                            strcmpi(D.mGoalcPars.GetValue(scan_sequence),'0') )
                        echo_pulse_sequence = 'SPIN';
                    end
                end
                vCurPar = 'EchoPulseSequence';  D.set_par(vCurPar, echo_pulse_sequence,'static');
                vCurPar = 'SteadyStatePulseSequence';  D.set_par(vCurPar, steady_state_pulse_seq,'static');
                
                if( strcmpi(acquisition_type, '3D'))
                    vCurPar = 'MRAcquisitionType';  D.set_par(vCurPar, '3D','static');
                else
                    vCurPar = 'MRAcquisitionType';  D.set_par(vCurPar, '2D','static');
                end
                
                vCurPar = 'MagneticFieldStrength';  D.set_par(vCurPar, D.mGoalcPars.GetValue(field_strength),'static');
                vCurPar = 'AcquisitionDuration';  D.set_par(vCurPar, D.mGoalcPars.GetValue(scan_duration),'static');
                
                % SpectrallySelectedExcitation
                if (strcmpi(D.mGoalcPars.GetValue(fat_saturation),'N') || ...
                        strcmpi(D.mGoalcPars.GetValue(fat_saturation),'0') )
                    SpectrallySelectedExcitation = 'NONE';
                else
                    SpectrallySelectedExcitation = 'FAT';
                end
                vCurPar = 'SpectrallySelectedExcitation';  D.set_par(vCurPar, SpectrallySelectedExcitation,'static');
                
                if( str2double(D.mGoalcPars.GetValue(phase_contrast_nr_dirs)) > 1 )
                    vCurPar = 'PhaseContrast';  D.set_par(vCurPar, 'YES','static');
                else
                    vCurPar = 'PhaseContrast';  D.set_par(vCurPar, 'NO','static');
                end
                if( str2double(D.mGoalcPars.GetValue(epi_factor)) > 1 )
                    vCurPar = 'EchoPlanarPulseSequence';  D.set_par(vCurPar, 'YES','static');
                else
                    vCurPar = 'EchoPlanarPulseSequence';  D.set_par(vCurPar, 'NO','static');
                end
                
                if (D.mGoalcPars.IsParameter('EX_ACQ_imaging_sequence'))
                    
                    switch (D.mGoalcPars.GetValue('VAL01_FFE_contrast_enh', [], 1))
                        case 0 %MGUFFE_CE_NO:
                            steady_state_pulse_seq = 'NONE';
                        case 1 % MGUFFE_CE_T1:
                            steady_state_pulse_seq = 'LONGITUDINAL';
                        case 2 % MGUFFE_CE_T2:
                            steady_state_pulse_seq = 'TIME_REVERSE';
                        case 3 % MGUFFE_CE_HYBRID:
                            steady_state_pulse_seq = 'FREE_PRECESSION';
                        otherwise
                            steady_state_pulse_seq = 'NONE';
                    end
                    
                    % Tagging
                    if (D.mGoalcPars.GetValue('VAL01_TAG_ctagging') ~= 0)
                        
                        switch (D.mGoalcPars.GetValue('VAL01_TAG_dimension', [], 1))
                            
                            case 0 % MGU_TAG_DIM_LINES:
                                tagging = 'LINE';
                            case 1 %MGU_TAG_DIM_GRID:
                                tagging = 'GRID';
                            case 2 %MGU_TAG_DIM_2XLINES:
                                tagging = 'GRID';
                            otherwise
                                tagging = 'NONE';
                        end
                    else
                        tagging = 'NONE';
                    end
                    
                    % MTC
                    switch (D.mGoalcPars.GetValue('VAL01_MTC_enable', [], 1))
                        case 0 % MPU_MTC_MODE_NO:
                            mtc = 'NONE';
                        case 1 %MPU_MTC_MODE_ON_RES:
                            mtc = 'ON_RESONANCE';
                        case 2 %MPU_MTC_MODE_OFF_RES:
                            mtc = 'OFF_RESONANCE';
                        case 3 %MPU_MTC_MODE_OFF_RES_MULTI_PULSE:
                            mtc = 'OFF_RESONANCE';
                        otherwise
                            mtc = 'NONE';
                    end
                    
                    % SpectrallySelectedExcitation
                    if (D.mGoalcPars.GetValue('VAL01_SPIR_mode') == 3) % MPUSPIR_MODE_PROSET)
                        SpectrallySelectedExcitation = 'NONE';
                    else
                        if (D.mGoalcPars.GetValue('VAL01_SPIR_suppr_tissue') == 0)
                            SpectrallySelectedExcitation = 'WATER';
                        else
                            SpectrallySelectedExcitation = 'FAT';
                        end
                    end
                    
                    % SpectrallySelectedSuppression
                    SpectrallySelectedSuppression = 'NONE';
                    water_suppress = false;
                    fat_suppress = false;
                    if (D.mGoalcPars.GetValue('VAL01_MTC_sos') == 1)
                        water_suppress = true;
                        fat_suppress = true;
                    end
                    
                    if (D.mGoalcPars.GetValue('VAL01_SPIR_enable') == 1)
                        if (D.mGoalcPars.GetValue('VAL01_SPIR_enable') == 1) %MPUSPIR_SUPPR_WATER
                            water_suppress = true;
                        else
                            fat_suppress = true;
                        end
                    end
                    
                    if (D.mGoalcPars.GetValue('VAL01_WATSUP_method') ~= 0)% MGUWATSUP_METHOD_NO
                        water_suppress = true;
                    end
                    
                    switch (D.mGoalcPars.GetValue('VAL01_WATSUP_in_base', [], 1))
                        case 2 % MPUWATSUP_INBASE_WATER:
                            water_suppress = true;
                        case 1 % MPUWATSUP_INBASE_FAT:
                            fat_suppress = true;
                        case 3 % MPUWATSUP_INBASE_WATER_FAT:
                            water_suppress = true;
                            fat_suppress = true;
                    end
                    
                    if (water_suppress && fat_suppress)
                        SpectrallySelectedSuppression = 'FAT_AND_WATER';
                    end
                    
                    if (D.mGoalcPars.GetValue('VAL01_WATSUP_method') ~= 0) %MGUWATSUP_METHOD_NO
                        water_suppress = true;
                    end
                    
                    switch (D.mGoalcPars.GetValue('VAL01_WATSUP_in_base', [], 1))
                        case 2 % MPUWATSUP_INBASE_WATER:
                            water_suppress = true;
                        case 1 % MPUWATSUP_INBASE_FAT:
                            fat_suppress = true;
                        case 3 %MPUWATSUP_INBASE_WATER_FAT:
                            water_suppress = true;
                            fat_suppress = true;
                    end
                    
                    if (water_suppress && ~fat_suppress)
                        SpectrallySelectedSuppression = 'WATER';
                    end
                    if (~water_suppress && fat_suppress)
                        SpectrallySelectedSuppression = 'FAT';
                    end
                    if (water_suppress && fat_suppress)
                        SpectrallySelectedSuppression = 'FAT_AND_WATER';
                    end
                    
                    % GeometryOfKSpaceTraversal
                    switch (D.mGoalcPars.GetValue('RC_k_space_traj_type', [], 1))
                        case 0
                            GeometryOfKSpaceTraversal = 'RECTILINEAR';
                        case 1
                            GeometryOfKSpaceTraversal = 'RADIAL';
                        case 2
                            GeometryOfKSpaceTraversal = 'SPIRAL';
                        case 3
                            GeometryOfKSpaceTraversal = 'PROPELLER';
                        otherwise
                            GeometryOfKSpaceTraversal = 'RECTILINEAR';
                    end
                    
                    % RectilinearPhaseEncodeReordering
                    if (D.mGoalcPars.GetValue('VAL01_ACQ_imaging_sequence') == 3) % MGUACQ_SEQ_FFE
                        outer_order = D.mGoalcPars.GetValue('VAL01_FFE_outer_profile_orders', 1, 1);
                    else
                        outer_order = D.mGoalcPars.GetValue('VAL01_RFE_outer_profile_orders', 1, 1);
                    end
                    
                    switch (outer_order)
                        case 2 % MGOBJORD_PROF_ORD_LINEAR:
                            RectilinearPhaseEncodeReordering = 'LINEAR';
                        case 3 % MGOBJORD_PROF_ORD_REV_LINEAR:
                            RectilinearPhaseEncodeReordering = 'REVERSE_LINEAR';
                        case 4 % MGOBJORD_PROF_ORD_SWEEP:
                            RectilinearPhaseEncodeReordering = 'SEGMENTED';
                        case 8 % MGOBJORD_PROF_ORD_SWEEP_2:
                            RectilinearPhaseEncodeReordering = 'SEGMENTED';
                        case 9 % MGOBJORD_PROF_ORD_SWEEP_4:
                            RectilinearPhaseEncodeReordering = 'SEGMENTED';
                        case 5 % MGOBJORD_PROF_ORD_OPT_CENTRIC:
                            RectilinearPhaseEncodeReordering = 'CENTRIC';
                        case 7 % MGOBJORD_PROF_ORD_OPT_MIN_CENTRIC:
                            RectilinearPhaseEncodeReordering = 'CENTRIC';
                        case 1 % MGOBJORD_PROF_ORD_LOW_HIGH:
                            RectilinearPhaseEncodeReordering = 'CENTRIC';
                        case 6 % MGOBJORD_PROF_ORD_OPT_REV_CENTRIC:
                            RectilinearPhaseEncodeReordering = 'REVERSE_CENTRIC';
                        case 0 % MGOBJORD_PROF_ORD_HIGH_LOW:
                            RectilinearPhaseEncodeReordering = 'REVERSE_CENTRIC';
                        case 10 % MGOBJORD_PROF_ORD_UNDEF:
                            RectilinearPhaseEncodeReordering = 'UNKNOWN';
                        case 12 % MGOBJORD_PROF_ORD_RADIAL:
                            RectilinearPhaseEncodeReordering = 'UNKNOWN';
                        case 13 % MGOBJORD_PROF_ORD_SPIRAL:
                            RectilinearPhaseEncodeReordering = 'UNKNOWN';
                        case 14 % MGOBJORD_PROF_ORD_VISTA:
                            RectilinearPhaseEncodeReordering = 'UNKNOWN';
                        otherwise
                            RectilinearPhaseEncodeReordering = 'UNKNOWN';
                    end
                    
                    % ParallelAcquisitionTechnique
                    ParallelAcquisitionTechnique = '';
                    if (D.mGoalcPars.GetValue('VAL01_GEO_sense_enable', [], 1) == 1)
                        ParallelAcquisitionTechnique = 'SENSE';
                    end
                    if (D.mGoalcPars.GetValue('VAL01_ACQ_kt_factor') > 1)
                        ParallelAcquisitionTechnique = 'KT';
                    end
                    if (D.mGoalcPars.IsParameter('VAL01_ACQ_MB') && (D.mGoalcPars.GetValue('VAL01_ACQ_MB', [], 1) == 1))
                        ParallelAcquisitionTechnique = 'MULTIBAND';
                    end
                    if (D.mGoalcPars.IsParameter('VAL01_ACQ_caipi') && (D.mGoalcPars.GetValue('VAL01_ACQ_caipi', [], 1) == 1))
                        ParallelAcquisitionTechnique = 'MULTIBAND';
                    end
                    
                    % CardiacSignalSource
                    CardiacSignalSource = '';
                    if (D.mGoalcPars.GetValue('VAL01_CARD_synchronisation', [], 1) ~= 0)
                        
                        switch (D.mGoalcPars.GetValue('EX_CARD_device', [], 1))
                            case 0 % MGG_CARD_DEV_ECG
                                CardiacSignalSource = 'ECG';
                            case 2 % MGG_CARD_DEV_INT
                                CardiacSignalSource = 'ECG';
                            case 3 % MGG_CARD_DEV_EXT
                                CardiacSignalSource = 'ECG';
                            case 1 % MGG_CARD_DEV_PPU
                                CardiacSignalSource = 'PP';
                        end
                    end
                    
                    % CoverageOfKSpace
                    if (D.mGoalcPars.GetValue('VAL01_ACQ_scan_mode', [], 1) == 1)
                        is_shutter = D.mGoalcPars.GetValue('VAL00_DEF_elliptical_k_space_shutter', [], 1) == 1;
                        if (is_shutter)
                            if (D.mGoalcPars.GetValue('VAL01_ACQ_scan_type', [], 1) == 1) % MGUACQ_SCT_SPECTRO
                                CoverageOfKSpace = 'ELLIPSOIDAL';
                            else
                                CoverageOfKSpace = 'ELLIPTICAL';
                            end
                            
                        else
                            CoverageOfKSpace = 'FULL';
                        end
                    else
                        CoverageOfKSpace = '';
                    end
                    
                    if (D.mGoalcPars.IsParameter('VAL01_DIFF_nr_directions'))
                        NrGradOris = D.mGoalcPars.GetValue('VAL01_DIFF_nr_directions');
                    elseif (D.mGoalcPars.IsParameter('VAL05_DIFF_measured_nr_oris'))
                        NrGradOris = D.mGoalcPars.GetValue('VAL05_DIFF_measured_nr_oris');
                    else
                        NrGradOris = 1;
                    end
                    
                    % PatientPosition
                    if( D.mGoalcPars.GetValue('EX_GEO_patient_position', [], 1) == 0)
                        patient_position = 'HF';
                    else
                        patient_position = 'FF';
                    end
                    if( D.mGoalcPars.GetValue('EX_GEO_patient_orientation', [], 1) == 0)
                        patient_position = [patient_position, 'S'];
                    elseif( D.mGoalcPars.GetValue('EX_GEO_patient_orientation', [], 1) == 1)
                        patient_position = [patient_position, 'P'];
                    elseif( D.mGoalcPars.GetValue('EX_GEO_patient_orientation', [], 1) == 2)
                        patient_position = [patient_position, 'L'];
                    else
                        patient_position = [patient_position, 'R'];
                    end
                    
                    vCurPar = 'PatientPosition'; D.set_par(vCurPar, patient_position ,'static');
                    
                    if( D.mGoalcPars.GetValue('VAL01_ACQ_scan_mode', [], 1) == 1 )
                        vCurPar = 'MRAcquisitionType';  D.set_par(vCurPar, '3D','static');
                    else
                        vCurPar = 'MRAcquisitionType';  D.set_par(vCurPar, '2D','static');
                    end
                    
                    vCurPar = 'MagneticFieldStrength';  D.set_par(vCurPar, D.mGoalcPars.GetValue('HW_main_magnetic_field_mT') / 1000,'static');
                    vCurPar = 'PercentPhaseFieldOfView';  D.set_par(vCurPar, D.mGoalcPars.GetValue('EX_GEO_rect_fov_perc'),'static');
                    vCurPar = 'AcquisitionDuration';  D.set_par(vCurPar, D.mGoalcPars.GetValue('AC_total_scan_time'),'static');
                    vCurPar = 'AcquisitionContrast';  D.set_par(vCurPar, D.mGoalcPars.GetValue('EX_FFE_ceffe'),'static');
                    
                    if( D.mGoalcPars.GetValue('VAL01_ACQ_fast_imaging_mode', [], 1) == 1 )
                        vCurPar = 'MultipleSpinEcho';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'MultipleSpinEcho';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    if( D.mGoalcPars.IsParameter('VAL01_ACQ_MB') && D.mGoalcPars.GetValue('VAL01_ACQ_MB', [], 1) == 1 )
                        vCurPar = 'MultiPlanarExcitation';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'MultiPlanarExcitation';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    if( D.mGoalcPars.GetValue('EX_PC_angio_mode', [], 1) == 2 )
                        vCurPar = 'PhaseContrast';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'PhaseContrast';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    if( D.mGoalcPars.GetValue('EX_PC_angio_mode', [], 1) == 1 )
                        vCurPar = 'TimeOfFlightContrast';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'TimeOfFlightContrast';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    if( D.mGoalcPars.GetValue('MPF_RFE_rf_spoil_mode', [], 1) == 1 && D.mGoalcPars.GetValue('VAL01_FFE_contrast_enh', [], 1) == 1 )
                        vCurPar = 'Spoiling';  D.set_par(vCurPar, 'RF','static');
                    else
                        vCurPar = 'Spoiling';  D.set_par(vCurPar, 'NONE','static');
                    end
                    
                    vCurPar = 'SteadyStatePulseSequence';  D.set_par(vCurPar, steady_state_pulse_seq,'static');
                    vCurPar = 'Tagging';  D.set_par(vCurPar, tagging,'static');
                    
                    if( D.mGoalcPars.GetValue('VAL01_ACQ_fast_imaging_mode', [], 1) == 3 || D.mGoalcPars.GetValue('VAL01_ACQ_fast_imaging_mode', [], 1) == 5 || D.mGoalcPars.GetValue('VAL01_ACQ_fast_imaging_mode', [], 1) == 4 )
                        vCurPar = 'EchoPlanarPulseSequence';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'EchoPlanarPulseSequence';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    vCurPar = 'MagnetizationTransfer';  D.set_par(vCurPar, mtc,'static');
                    
                    if( D.mGoalcPars.GetValue('VAL01_T2PREP_enable', [], 1) == 1 )
                        vCurPar = 'T2Preparation';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'T2Preparation';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    if( D.mGoalcPars.GetValue('VAL01_TFEPP_type', [], 1) == 3 )
                        vCurPar = 'BloodSignalNulling';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'BloodSignalNulling';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    if( D.mGoalcPars.GetValue('VAL01_ACQ_tfe_excitations', [], 1) == 1 && D.mGoalcPars.GetValue('VAL01_TFEPP_type', [], 1) == 2 )
                        vCurPar = 'SaturationRecovery';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'SaturationRecovery';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    vCurPar = 'SpectrallySelectedExcitation';  D.set_par(vCurPar, SpectrallySelectedExcitation,'static');
                    vCurPar = 'SpectrallySelectedSuppression';  D.set_par(vCurPar, SpectrallySelectedSuppression,'static');
                    vCurPar = 'GeometryOfKSpaceTraversal';  D.set_par(vCurPar, GeometryOfKSpaceTraversal,'static');
                    
                    if( D.mGoalcPars.GetValue('RC_quad_receive_coils', [], 1) ~= 0 )
                        vCurPar = 'QuadratureReceiveCoil';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'QuadratureReceiveCoil';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    vCurPar = 'RectilinearPhaseEncodeReordering';  D.set_par(vCurPar, RectilinearPhaseEncodeReordering,'static');
                    
                    if( D.mGoalcPars.GetValue('RC_receive_coil_types', [], 1) == 3 )
                        vCurPar = 'ReceiveCoilType';  D.set_par(vCurPar, 'MULTICOIL','static');
                    else
                        if( D.mGoalcPars.GetValue('RC_receive_coil_types', [], 1) == 4 )
                            vCurPar = 'ReceiveCoilType';  D.set_par(vCurPar, 'UNKNOWN','static');
                        else
                            vCurPar = 'ReceiveCoilType';  D.set_par(vCurPar, strrep(D.mGoalcPars.GetValue('RC_receive_coil_types', 1), 'ARGRC_COIL_TYPE_', ''),'static');
                        end
                    end
                    
                    if( D.mGoalcPars.GetValue('VAL01_SPY_decouple', [], 1) ~= 0 )
                        vCurPar = 'Decoupling';  D.set_par(vCurPar, 'YES','static');
                    else
                        vCurPar = 'Decoupling';  D.set_par(vCurPar, 'NO','static');
                    end
                    
                    vCurPar = 'ParallelAcquisitionTechnique';  D.set_par(vCurPar, ParallelAcquisitionTechnique,'static');
                    vCurPar = 'CardiacSignalSource';  D.set_par(vCurPar, CardiacSignalSource,'static');
                    vCurPar = 'CoverageOfKSpace';  D.set_par(vCurPar, CoverageOfKSpace,'static');
                    vCurPar = {'0x2001', '0x1009'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('VW_prepulse_delay'),'static');
                    vCurPar = {'0x2005', '0x1415'};  D.set_par(vCurPar, NrGradOris,'static');
                    
                    vCurPar = {'0x2005', '0x1054'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_ap_angulations', 1,'static'));
                    vCurPar = {'0x2005', '0x1055'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_fh_angulations', 1,'static'));
                    vCurPar = {'0x2005', '0x1056'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_rl_angulations', 1,'static'));
                    vCurPar = {'0x2005', '0x1057'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_ap_fovs', 1,'static'));
                    vCurPar = {'0x2005', '0x1058'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_fh_fovs', 1,'static'));
                    vCurPar = {'0x2005', '0x1059'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_rl_fovs', 1,'static'));
                    vCurPar = {'0x2005', '0x105a'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_ap_offcentres', 1,'static'));
                    vCurPar = {'0x2005', '0x105b'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_fh_offcentres', 1,'static'));
                    vCurPar = {'0x2005', '0x105c'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_volume_rl_offcentres', 1,'static'));
                    vCurPar = {'0x2005', '0x105e'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('PS_slab_orientations', 1,'static'));
                    vCurPar = {'0x2005', '0x105d'};  D.set_par(vCurPar, upper(D.mGoalcPars.GetValue('PS_volume_geom_types', 1,'static')));
                    vCurPar = {'0x2001', '0x1001'};  D.set_par(vCurPar, D.mGoalcPars.GetValue('RC_chemical_shifts', 1,'static'));
                end
            catch
                warning('Could not get all DICOM tags. Please make sure you have the most recent version of the ReconFrame patch');
            end
        end
        function dicom_get_static_values_from_parfile(D, Parfile)
            date_time = DICOMExporter.format_date_time(Parfile.ExaminationDateTime);
            if( ~isempty(date_time ))
                if( ~isempty(date_time{1} ))
                    
                    vCurPar = 'StudyDate'; D.set_par(vCurPar, date_time{1},'static');
                    vCurPar = 'StudyDate'; D.set_par(vCurPar, date_time{1},'static');
                    vCurPar = 'SeriesDate'; D.set_par(vCurPar, date_time{1},'static');
                    vCurPar = 'AcquisitionDate'; D.set_par(vCurPar, date_time{1},'static');
                    vCurPar = 'ContentDate'; D.set_par(vCurPar, date_time{1},'static');
                    vCurPar = 'AcquisitionDateTime'; D.set_par(vCurPar, date_time{1},'static');
                end
                if( ~isempty(date_time{2} ))
                    vCurPar = 'StudyTime'; D.set_par(vCurPar, date_time{2},'static');
                    vCurPar = 'AcquisitionTime'; D.set_par(vCurPar, date_time{2},'static');
                    vCurPar = 'ContentTime'; D.set_par(vCurPar, date_time{2},'static');
                end
            end
            
            vCurPar = 'Modality'; D.set_par(vCurPar, 'MR','static');
            vCurPar = 'ProtocolName'; D.set_par(vCurPar, Parfile.ProtocolName,'static');
            vCurPar = {'0x2005', '0x1035'}; D.set_par(vCurPar, Parfile.SeriesDataType,'static');
            vCurPar = 'AcquisitionNumber'; D.set_par(vCurPar, Parfile.AcquisitionNr,'static');
            vCurPar = {'0x2001', '0x101d'}; D.set_par(vCurPar, Parfile.ReconstructionNr,'static');
            vCurPar = 'AcquisitionDuration'; D.set_par(vCurPar, Parfile.ScanDuration,'static');
            vCurPar = {'0x2001', '0x1017'}; D.set_par(vCurPar, Parfile.MaxNumberOfCardiacPhases,'static');
            vCurPar = {'0x2001', '0x1014'}; D.set_par(vCurPar, Parfile.MaxNumberOfEchoes,'static');
            vCurPar = {'0x2001', '0x1015'}; D.set_par(vCurPar, Parfile.MaxNumberOfSlicesLocations,'static');
            vCurPar = {'0x2001', '0x1081'}; D.set_par(vCurPar, Parfile.MaxNumberOfDynamics,'static');
            vCurPar = 'NumberOfTemporalPositions'; D.set_par(vCurPar, Parfile.MaxNumberOfDynamics,'static');
            vCurPar = 'PatientPosition'; D.set_par(vCurPar, strrep(Parfile.PatientPosition, '.', ''),'static');
            vCurPar = {'0x2005', '0x107b'}; D.set_par(vCurPar, Parfile.PreparationDirection,'static');
            vCurPar = {'0x2001', '0x1020'}; D.set_par(vCurPar, Parfile.Technique,'static');
            vCurPar = 'PulseSequenceName'; D.set_par(vCurPar, Parfile.Technique,'static');
            
            acquisition_matrix = [0; Parfile.ScanResolution(1); Parfile.ScanResolution(2); 0];
            vCurPar = 'AcquisitionMatrix'; D.set_par(vCurPar, acquisition_matrix,'static');
            vCurPar = 'MRAcquisitionFrequencyEncodingSteps'; D.set_par(vCurPar, Parfile.ScanResolution(1),'static');
            vCurPar = {'0x2005', '0x101D'}; D.set_par(vCurPar, Parfile.ScanResolution(1),'static');
            vCurPar = {'0x2005', '0x106f'}; D.set_par(vCurPar, Parfile.ScanMode,'static');
            vCurPar = 'RepetitionTime'; D.set_par(vCurPar, Parfile.RepetitionTime,'static');
            
            repetition_times = [Parfile.RepetitionTime; 0];
            vCurPar = {'0x2005', '0x1030'}; D.set_par(vCurPar, repetition_times,'static');
            vCurPar = {'0x2005', '0x1074'}; D.set_par(vCurPar, Parfile.FOV(1),'static');
            vCurPar = {'0x2005', '0x1075'}; D.set_par(vCurPar, Parfile.FOV(2),'static');
            vCurPar = {'0x2005', '0x1076'}; D.set_par(vCurPar, Parfile.FOV(3),'static');
            vCurPar = {'0x2001', '0x1022'}; D.set_par(vCurPar, Parfile.WaterFatShift,'static');
            vCurPar = {'0x2005', '0x1071'}; D.set_par(vCurPar, Parfile.AngulationMidslice(1),'static');
            vCurPar = {'0x2005', '0x1072'}; D.set_par(vCurPar, Parfile.AngulationMidslice(2),'static');
            vCurPar = {'0x2005', '0x1073'}; D.set_par(vCurPar, Parfile.AngulationMidslice(3),'static');
            vCurPar = {'0x2005', '0x1078'}; D.set_par(vCurPar, Parfile.OffCentreMidslice(1),'static');
            vCurPar = {'0x2005', '0x1079'}; D.set_par(vCurPar, Parfile.OffCentreMidslice(2),'static');
            vCurPar = {'0x2005', '0x107a'}; D.set_par(vCurPar, Parfile.OffCentreMidslice(3),'static');
            
            YES_NO = {'N', 'Y'};
            NONE_SLAB = {'NONE', 'SLAB'};
            vCurPar = {'0x2005', '0x1016'}; D.set_par(vCurPar, YES_NO{Parfile.FlowCompensation+1},'static');
            vCurPar = 'SpatialPresaturation'; D.set_par(vCurPar, NONE_SLAB{Parfile.Presaturation+1},'static');
            
            if( Parfile.ImageInformation.SliceOrientation(1) == 1 ) && ( strcmpi(strtrim(Parfile.PreparationDirection) , 'RL'))  % TRA, phase-enc in RL
                venc = [Parfile.PhaseEncodingVelocity(1);Parfile.PhaseEncodingVelocity(2);Parfile.PhaseEncodingVelocity(3)];
            elseif( Parfile.ImageInformation.SliceOrientation(1) == 1 ) && ( strcmpi(strtrim(Parfile.PreparationDirection) , 'AP'))   % TRA, phase-enc in AP  
                venc = [Parfile.PhaseEncodingVelocity(2);Parfile.PhaseEncodingVelocity(1);Parfile.PhaseEncodingVelocity(3)];        
            elseif( Parfile.ImageInformation.SliceOrientation(1) == 2 )  && ( strcmpi(strtrim(Parfile.PreparationDirection) , 'AP'))  % SAG, phase-enc in AP  
                venc = [Parfile.PhaseEncodingVelocity(1);Parfile.PhaseEncodingVelocity(2);Parfile.PhaseEncodingVelocity(3)];
            elseif( Parfile.ImageInformation.SliceOrientation(1) == 2 )  && ( strcmpi(strtrim(Parfile.PreparationDirection) , 'FH'))  % SAG, phase-enc in FH  
                venc = [Parfile.PhaseEncodingVelocity(2);Parfile.PhaseEncodingVelocity(1);Parfile.PhaseEncodingVelocity(3)];
            elseif( Parfile.ImageInformation.SliceOrientation(1) == 3 ) && ( strcmpi(strtrim(Parfile.PreparationDirection) , 'FH'))   % COR, phase-enc in FH
                venc = [Parfile.PhaseEncodingVelocity(2);Parfile.PhaseEncodingVelocity(1);Parfile.PhaseEncodingVelocity(3)];
            elseif( Parfile.ImageInformation.SliceOrientation(1) == 3 ) && ( strcmpi(strtrim(Parfile.PreparationDirection) , 'RL'))   % COR, phase-enc in RL
                venc = [Parfile.PhaseEncodingVelocity(1);Parfile.PhaseEncodingVelocity(2);Parfile.PhaseEncodingVelocity(3)];
            end
            vCurPar = {'0x2001', '0x101a'}; D.set_par(vCurPar, venc,'static');
            
            NONE_RESONANCE = {'NONE', 'ON_RESONANCE'};
            NONE_WATER = {'NONE', 'WATER'};
            vCurPar = 'MagnetizationTransfer'; D.set_par(vCurPar, NONE_RESONANCE{Parfile.MTC+1},'static');
            vCurPar = 'SpectrallySelectedExcitation'; D.set_par(vCurPar, NONE_WATER{Parfile.SPIR+1},'static');
            vCurPar = {'0x2001', '0x1013'}; D.set_par(vCurPar, Parfile.EPIFactor,'static');
            vCurPar = {'0x2001', '0x1012'}; D.set_par(vCurPar, YES_NO{Parfile.DynamicScan+1},'static');
            vCurPar = {'0x2005', '0x1014'}; D.set_par(vCurPar, YES_NO{Parfile.Diffusion+1},'static');
            vCurPar = {'0x2001', '0x1011'}; D.set_par(vCurPar, Parfile.DiffusionEchoTime,'static');
            vCurPar = {'0x2005', '0x1428'}; D.set_par(vCurPar, Parfile.NumberOfLabelTypes,'static');
            
            frame_of_reference_uid = dicomuid;
            study_uid = dicomuid;
            time_format = 'HHMMSS.FFF';
            date_format = 'yyyymmdd';
            cur_date = datestr(now,date_format);
            cur_time = datestr(now,time_format);
            
            vCurPar = 'FrameOfReferenceUID'; D.set_par(vCurPar, frame_of_reference_uid,'static');
            vCurPar = 'SOPClassUID'; D.set_par(vCurPar, '1.2.840.10008.5.1.4.1.1.4','static');
            vCurPar = 'MediaStorageSOPClassUID'; D.set_par(vCurPar, '1.2.840.10008.5.1.4.1.1.4','static');
            vCurPar = 'StudyInstanceUID'; D.set_par(vCurPar, study_uid,'static');
            vCurPar = 'InstanceCreationDate'; D.set_par(vCurPar, cur_date,'static');
            vCurPar = 'InstanceCreationTime'; D.set_par(vCurPar, cur_time,'static');
            vCurPar = 'SpecificCharacterSet'; D.set_par(vCurPar, 'ISO_IR 100','static');
            
            vCurPar = 'PatientID'; D.set_par(vCurPar, D.mPatientID,'static');
            
            % Calculate the Study ID. The Study ID is a short string and
            % cannot be longer than 16 chars. Therefore take the last 16
            % characters from the StudyUID
            vCurPar = 'StudyID'; D.set_par(vCurPar, study_uid(end-15:end),'static');
                        
            vCurPar = 'SeriesNumber'; D.set_par(vCurPar, Parfile.AcquisitionNr,'static');
            vCurPar = 'PatientBirthDate'; D.set_par(vCurPar, '20000101','static');
            vCurPar = 'PatientSex'; D.set_par(vCurPar, 'O','static');
            vCurPar = 'PatientName'; D.set_par(vCurPar, Parfile.PatientName,'static');
            vCurPar = 'InstitutionName'; D.set_par(vCurPar, 'NA','static');
            vCurPar = 'InstitutionAddress'; D.set_par(vCurPar, 'NA','static');
            vCurPar = 'PerformingPhysicianName'; D.set_par(vCurPar, 'NA','static');
            vCurPar = 'Manufacturer'; D.set_par(vCurPar, 'Philips Healthcare','static');
            
            vCurPar = 'CodingSchemeDesignator'; D.set_par(vCurPar, 'DCM','static');
            vCurPar = 'SamplesPerPixel'; D.set_par(vCurPar, 1,'static');
            vCurPar = 'VolumetricProperties'; D.set_par(vCurPar, 'VOLUME','static');
            vCurPar = 'VolumeBasedCalculationTechnique'; D.set_par(vCurPar, 'NONE','static');
            vCurPar = 'TimeDomainFiltering'; D.set_par(vCurPar, '','static');
            vCurPar = 'MetaboliteMapDescription'; D.set_par(vCurPar, 'WATER','static');
            vCurPar = 'KSpaceFiltering'; D.set_par(vCurPar, 'HAMMING','static');
            vCurPar = 'LUTLabel'; D.set_par(vCurPar, 'Philips','static');
            vCurPar = {'0x2005', '0x1035'}; D.set_par(vCurPar, 'PIXEL','static');
            vCurPar = {'0x2001', '0x101d'}; D.set_par(vCurPar, -1,'static');
            
            % still empty (no CRecon parameter found)
            vCurPar = 'ReferringPhysicianName'; D.set_par(vCurPar, '','static');
            vCurPar = 'AdditionalPatientHistory'; D.set_par(vCurPar, '','static');
            vCurPar = 'SequenceVariant'; D.set_par(vCurPar, '','static');
            vCurPar = 'ScanOptions'; D.set_par(vCurPar, '','static');
            vCurPar = 'MultiCoilElementName'; D.set_par(vCurPar, '','static');
            vCurPar = 'DecoupledNucleus'; D.set_par(vCurPar, '','static');
            
            % Harcoded values (TODO: find the proper values for them)
            vCurPar = {'0x2005', '0x1061'}; D.set_par(vCurPar, 'NO','static');	% Image Prepulse
            vCurPar = 'DataPointRows'; D.set_par(vCurPar, 1,'static');
            vCurPar = 'DataPointColumns'; D.set_par(vCurPar, 0,'static');
            
            % Add the PrivateCreator's
            vCurPar = 'PrivateCreator110'; D.set_par(vCurPar, 'Philips Imaging DD 001','static');
            vCurPar = 'PrivateCreator111'; D.set_par(vCurPar, 'Philips Imaging DD 002','static');
            
            vCurPar = 'PrivateCreator510'; D.set_par(vCurPar, 'Philips MR Imaging DD 001','static');
            vCurPar = 'PrivateCreator511'; D.set_par(vCurPar, 'Philips MR Imaging DD 002','static');
            vCurPar = 'PrivateCreator512'; D.set_par(vCurPar, 'Philips MR Imaging DD 003','static');
            vCurPar = 'PrivateCreator513'; D.set_par(vCurPar, 'Philips MR Imaging DD 004','static');
            vCurPar = 'PrivateCreator514'; D.set_par(vCurPar, 'Philips MR Imaging DD 005','static');
            
            vCurPar = 'TransferSyntaxUID'; D.set_par(vCurPar, '1.2.840.10008.1.2.1','static');
        end
        function dicom_get_database_values(D)
            if( ~isempty(D.mGoalcPars) )
                if (strcmpi(D.mRelease(1:4), 'REL5'))
                    % Automatically generated by Matlab script format_dicom_tags.m
                    vCurPar = {'0x2005', '0x143A'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_DATA_DICTIONARY_CONTENTS_VERSION') ,'static');
                    vCurPar = {'0x2005', '0x1397'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_ANATOMIC_REG_CODE_VALUE') ,'static');
                    vCurPar = {'0x2005', '0x1249'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_MR_NR_OF_SERIES_OPERATORS_NAME') ,'static');
                    vCurPar = {'0x2005', '0x1245'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_MR_NR_OF_SOFTWARE_VERSION') ,'static');
                    vCurPar = {'0x2005', '0x1037'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_MR_SERIES_IS_SPECTRO') ,'static');
                    vCurPar = {'0x2001', '0x10CC'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_SERIES_DERIVATION_DESCRIPTION') ,'static');
                    vCurPar = {'0x2001', '0x107B'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_MR_SERIES_ACQUISITION_NUMBER') ,'static');
                    vCurPar = {'0x2001', '0x106E'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_SERIES_TYPE') ,'static');
                    vCurPar = {'0x2001', '0x1062'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_SERIES_COMMITTED') ,'static');
                    vCurPar = {'0x2001', '0x1061'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_SERIES_TRANSMITTED') ,'static');
                    vCurPar = {'0x1001', '0x100C'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_REMOTE_ACTIVITY_STATUS') ,'static');
                    vCurPar = {'0x1001', '0x100B'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_LOCAL_MEDIA_WRITE_STATUS') ,'static');
                    vCurPar = {'0x1001', '0x100A'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_STORAGE_COMMIT_STATUS') ,'static');
                    vCurPar = {'0x1001', '0x1009'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_ARCHIVE_TRANSFER_STATUS') ,'static');
                    vCurPar = {'0x1001', '0x1008'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_EXPORT_STATUS') ,'static');
                    vCurPar = {'0x1001', '0x1007'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_PIIM_PRINT_STATUS') ,'static');
                    vCurPar = {'0x0020', '0x0060'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_DICOM_LATERALITY') ,'static');
                    vCurPar = {'0x0020', '0x0052'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_DICOM_FRAME_OF_REFERENCE_UID') ,'static');
                    vCurPar = {'0x0008', '0x0021'}; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_DICOM_SERIES_DATE', true) ,'static');
                    vCurPar = {'0x0018', '0x9089'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_DIFFUSION_GRADIENT_ORIENTATION') ,'static');
                    vCurPar = {'0x0018', '0x9087'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_DIFFUSION_BVALUE') ,'static');
                    vCurPar = {'0x0018', '0x9080'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_METABOLITE_MAP_DESCRIPTION') ,'static');
                    vCurPar = {'0x0018', '0x9079'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_INVERSION_TIMES') ,'static');
                    vCurPar = {'0x0018', '0x9075'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_DIFFUSION_DIRECTIONALITY') ,'static');
                    vCurPar = {'0x0018', '0x9064'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_K_SPACE_FILTERING') ,'static');
                    vCurPar = {'0x0018', '0x9047'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_MULTI_COIL_ELEMENT_NAME') ,'static');
                    vCurPar = {'0x0018', '0x9044'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_QUADRATURE_RECEIVE_COIL') ,'static');
                    vCurPar = {'0x0018', '0x9043'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_RECEIVE_COIL_TYPE') ,'static');
                    vCurPar = {'0x0018', '0x1314'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_FLIP_ANGLE') ,'static');
                    vCurPar = {'0x0018', '0x1312'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_PHASE_ENCODING_DIRECTION') ,'static');
                    vCurPar = {'0x0018', '0x1310'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_ACQUISITION_MATRIX') ,'static');
                    vCurPar = {'0x0018', '0x1250'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_RECEIVING_COIL') ,'static');
                    vCurPar = {'0x0018', '0x1094'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_TRIGGER_WINDOW') ,'static');
                    vCurPar = {'0x0018', '0x0080'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_REPETITION_TIME') ,'static');
                    vCurPar = {'0x0018', '0x0050'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_SLICE_THICKNESS') ,'static');
                    vCurPar = {'0x0018', '0x0023'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_MR_ACQUISITION_TYPE') ,'static');
                    vCurPar = {'0x0018', '0x0022'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_SCAN_OPTIONS') ,'static');
                    vCurPar = {'0x0018', '0x0021'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_SEQUENCE_VARIANT') ,'static');
                    vCurPar = {'0x0018', '0x0020'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_SCANNING_SEQUENCE') ,'static');
                    vCurPar = {'0x0018', '0x0010'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_CONTRAST_BOLUS_AGENT') ,'static');
                    vCurPar = {'0x0008', '0x9209'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_ACQUISITION_CONTRAST') ,'static');
                    vCurPar = {'0x0008', '0x0033'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_CONTENT_TIME') ,'static');
                    vCurPar = {'0x0008', '0x002A'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_ACQUISITION_DATETIME', true) ,'static');
                    vCurPar = {'0x0008', '0x0023'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_CONTENT_DATE', true) ,'static');
                    vCurPar = {'0x0008', '0x0022'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_ACQUISITION_DATE', true) ,'static');
                    vCurPar = {'0x0008', '0x0014'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_INSTANCE_CREATOR_UID') ,'static');
                    vCurPar = {'0x0008', '0x0008'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_IMAGE_TYPE') ,'static');
                    vCurPar = {'0x0028', '0x0011'}; D.set_par(vCurPar, D.get_db_value('RFR_FOREIGN_DICOM_COLUMNS') ,'static');
                    vCurPar = {'0x0028', '0x0010'}; D.set_par(vCurPar, D.get_db_value('RFR_FOREIGN_DICOM_ROWS') ,'static');
                    vCurPar = {'0x0018', '0x9044'}; D.set_par(vCurPar, D.get_db_value('RFR_FOREIGN_DICOM_QUADRATURE_RECEIVE_COIL') ,'static');
                    vCurPar = {'0x0018', '0x1094'}; D.set_par(vCurPar, D.get_db_value('RFR_FOREIGN_DICOM_TRIGGER_WINDOW') ,'static');
                    vCurPar = {'0x0018', '0x1088'}; D.set_par(vCurPar, D.get_db_value('RFR_FOREIGN_DICOM_HEART_RATE') ,'static');
                    vCurPar = {'0x0018', '0x1084'}; D.set_par(vCurPar, D.get_db_value('RFR_FOREIGN_DICOM_INTERVALS_REJECTED') ,'static');
                    vCurPar = {'0x0018', '0x1083'}; D.set_par(vCurPar, D.get_db_value('RFR_FOREIGN_DICOM_INTERVALS_ACQUIRED') ,'static');
                    vCurPar = {'0x0018', '0x1082'}; D.set_par(vCurPar, D.get_db_value('RFR_FOREIGN_DICOM_HIGH_RR_VALUE') ,'static');
                    vCurPar = {'0x0018', '0x1081'}; D.set_par(vCurPar, D.get_db_value('RFR_FOREIGN_DICOM_LOW_RR_VALUE') ,'static');
                    vCurPar = {'0x0018', '0x0094'}; D.set_par(vCurPar, D.get_db_value('RFR_FOREIGN_DICOM_PERCENT_PHASE_FIELD_OF_VIEW') ,'static');
                    vCurPar = {'0x0018', '0x0093'}; D.set_par(vCurPar, D.get_db_value('RFR_FOREIGN_DICOM_PERCENT_SAMPLING') ,'static');
                    vCurPar = {'0x0018', '0x0091'}; D.set_par(vCurPar, D.get_db_value('RFR_FOREIGN_DICOM_ECHO_TRAIN_LENGTH') ,'static');
                    vCurPar = {'0x0018', '0x0089'}; D.set_par(vCurPar, D.get_db_value('RFR_FOREIGN_DICOM_NUMBER_OF_PHASE_ENCODING_STEPS') ,'static');
                    vCurPar = {'0x0018', '0x0082'}; D.set_par(vCurPar, D.get_db_value('RFR_FOREIGN_DICOM_INVERSION_TIME') ,'static');
                    vCurPar = {'0x0018', '0x0050'}; D.set_par(vCurPar, D.get_db_value('RFR_FOREIGN_DICOM_SLICE_THICKNESS') ,'static');
                    vCurPar = {'0x2005', '0x1429'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_IMAGE_LABEL_TYPE') ,'static');
                    vCurPar = {'0x2005', '0x1413'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_IMAGE_GRADIENT_ORIENTATION_NUMBER') ,'static');
                    vCurPar = {'0x2005', '0x1412'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_IMAGE_DIFF_B_VALUE_NUMBER') ,'static');
                    vCurPar = {'0x2005', '0x10B2'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_IMAGE_DIFFUSION_FH') ,'static');
                    vCurPar = {'0x2005', '0x10B1'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_IMAGE_DIFFUSION_AP') ,'static');
                    vCurPar = {'0x2005', '0x10B0'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_IMAGE_DIFFUSION_RL') ,'static');
                    vCurPar = {'0x2005', '0x10A8'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_PRIVATE_INVERSION_TIME') ,'static');
                    vCurPar = {'0x2005', '0x10A1'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_SYNCRA_SCAN_TYPE') ,'static');
                    vCurPar = {'0x2005', '0x10A0'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_DYNAMIC_SCAN_BEGIN_TIME') ,'static');
                    vCurPar = {'0x2005', '0x106E'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_IMAGE_SCANNING_SEQUENCE_PRIVATE') ,'static');
                    vCurPar = {'0x2005', '0x1063'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_M_RF_MRI_STATUS_INDICATION') ,'static');
                    vCurPar = {'0x2005', '0x1061'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_IMAGE_PREPULSE_TYPE') ,'static');
                    vCurPar = {'0x2005', '0x1011'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_IMAGE_TYPE_MR') ,'static');
                    vCurPar = {'0x2005', '0x100E'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_SCALE_SLOPE') ,'static');
                    vCurPar = {'0x2005', '0x100D'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_SCALE_INTERCEPT') ,'static');
                    vCurPar = {'0x2005', '0x100C'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_MIN_FP') ,'static');
                    vCurPar = {'0x2005', '0x100B'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_MAX_FP') ,'static');
                    vCurPar = {'0x2001', '0x1003'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_IMAGE_DIFFUSION_B_FACTOR') ,'static');
                    vCurPar = {'0x2001', '0x1002'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_IMAGE_CHEMICAL_SHIFT_NUMBER') ,'static');
                    vCurPar = {'0x2001', '0x1001'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_PIIM_MR_IMAGE_CHEMICAL_SHIFT') ,'static');
%                     vCurPar = {'0x0028', '0x1054'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_RESCALE_TYPE') ,'static');
%                     vCurPar = {'0x0028', '0x1053'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_RESCALE_SLOPE') ,'static');
%                     vCurPar = {'0x0028', '0x1052'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_RESCALE_INTERCEPT') ,'static');
%                     vCurPar = {'0x2005', '0x140B'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_RESCALE_TYPE') ,'static');       % LMG edit: add dicom tag for circle
%                     vCurPar = {'0x2005', '0x140A'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_RESCALE_SLOPE') ,'static');      % LMG edit: add dicom tag for circle
%                     vCurPar = {'0x2005', '0x1409'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_RESCALE_INTERCEPT') ,'static');  % LMG edit: add dicom tag for circle
                    vCurPar = {'0x0028', '0x1051'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_WINDOW_WIDTH') ,'static');
                    vCurPar = {'0x0028', '0x1050'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_WINDOW_CENTER') ,'static');
                    vCurPar = {'0x0028', '0x0103'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_PIXEL_REPRESENTATION') ,'static');
                    vCurPar = {'0x0028', '0x0102'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_HIGH_BIT') ,'static');
                    vCurPar = {'0x0028', '0x0101'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_BITS_STORED') ,'static');
                    vCurPar = {'0x0028', '0x0100'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_BITS_ALLOCATED') ,'static');
                    vCurPar = {'0x0028', '0x0034'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_PIXEL_ASPECT_RATIO') ,'static');
                    vCurPar = {'0x0028', '0x0030'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_PIXEL_SPACING') ,'static');
                    vCurPar = {'0x0028', '0x0011'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_COLUMNS') ,'static');
                    vCurPar = {'0x0028', '0x0010'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_ROWS') ,'static');
                    vCurPar = {'0x0028', '0x0006'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_PLANAR_CONFIGURATION') ,'static');
                    vCurPar = {'0x0028', '0x0004'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_PHOTOMETRIC_INTERPRETATION') ,'static');
                    vCurPar = {'0x0028', '0x0002'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_SAMPLES_PER_PIXEL') ,'static');
                    vCurPar = {'0x0020', '0x1041'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_SLICE_LOCATION') ,'static');
                    vCurPar = {'0x0020', '0x0105'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_NUMBER_OF_TEMPORAL_POSITIONS') ,'static');
                    vCurPar = {'0x0020', '0x0100'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_TEMPORAL_POSITION_IDENTIFIER') ,'static');
                    vCurPar = {'0x0018', '0x9603'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_DIFFUSSION_B_VALUE_XY') ,'static');
                    vCurPar = {'0x0018', '0x9602'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_DIFFUSSION_B_VALUE_XX') ,'static');
                    vCurPar = {'0x0018', '0x9232'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_MR_ACQUISITION_PHASE_ENCODING_STEPS_IOU_OF_PLANE') ,'static');
                    vCurPar = {'0x0018', '0x9147'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_DIFFUSION_ANISOTROPY_TYPE') ,'static');
                    vCurPar = {'0x0018', '0x9089'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_DIFFUSION_GRADIENT_ORIENTATION') ,'static');
                    vCurPar = {'0x0018', '0x9087'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_DIFFUSION_BVALUE') ,'static');
                    vCurPar = {'0x0018', '0x9080'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_METABOLITE_MAP_DESCRIPTION') ,'static');
                    vCurPar = {'0x0018', '0x9079'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_INVERSION_TIMES') ,'static');
                    vCurPar = {'0x0018', '0x9075'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_DIFFUSION_DIRECTIONALITY') ,'static');
                    vCurPar = {'0x0018', '0x9064'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_K_SPACE_FILTERING') ,'static');
                    vCurPar = {'0x0018', '0x9047'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_MULTI_COIL_ELEMENT_NAME') ,'static');
                    vCurPar = {'0x0018', '0x9044'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_QUADRATURE_RECEIVE_COIL') ,'static');
                    vCurPar = {'0x0018', '0x9043'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_RECEIVE_COIL_TYPE') ,'static');
                    vCurPar = {'0x0018', '0x1314'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_FLIP_ANGLE') ,'static');
                    vCurPar = {'0x0018', '0x1312'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_PHASE_ENCODING_DIRECTION') ,'static');
                    vCurPar = {'0x0018', '0x1310'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_ACQUISITION_MATRIX') ,'static');
                    vCurPar = {'0x0018', '0x1250'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_RECEIVING_COIL') ,'static');
                    vCurPar = {'0x0018', '0x1094'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_TRIGGER_WINDOW') ,'static');
                    vCurPar = {'0x0018', '0x0080'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_REPETITION_TIME') ,'static');
                    vCurPar = {'0x0018', '0x0050'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_SLICE_THICKNESS') ,'static');
                    vCurPar = {'0x0018', '0x0023'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_MR_ACQUISITION_TYPE') ,'static');
                    vCurPar = {'0x0018', '0x0022'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_SCAN_OPTIONS') ,'static');
                    vCurPar = {'0x0018', '0x0021'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_SEQUENCE_VARIANT') ,'static');
                    vCurPar = {'0x0018', '0x0020'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_SCANNING_SEQUENCE') ,'static');
                    vCurPar = {'0x0018', '0x0010'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_CONTRAST_BOLUS_AGENT') ,'static');
                    vCurPar = {'0x0008', '0x9209'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_ACQUISITION_CONTRAST') ,'static');
                    vCurPar = {'0x0008', '0x0033'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_CONTENT_TIME') ,'static');
                    vCurPar = {'0x0008', '0x002A'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_ACQUISITION_DATETIME', true) ,'static');
                    vCurPar = {'0x0008', '0x0023'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_CONTENT_DATE', true) ,'static');
                    vCurPar = {'0x0008', '0x0022'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_ACQUISITION_DATE', true) ,'static');
                    vCurPar = {'0x0008', '0x0014'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_INSTANCE_CREATOR_UID') ,'static');
                    vCurPar = {'0x0008', '0x0008'}; D.set_par(vCurPar, D.get_db_value('RFR_IMAGE_DICOM_IMAGE_TYPE') ,'static');
                    vCurPar = {'0x2005', '0x1147'}; D.set_par(vCurPar, D.get_db_value('RFR_MRBLOBDATA_PIIM_MR_BLOB_FLAG') ,'static');
                    vCurPar = {'0x2005', '0x1144'}; D.set_par(vCurPar, D.get_db_value('RFR_MRBLOBDATA_PIIM_MR_BLOB_DATA') ,'static');
                    vCurPar = {'0x2005', '0x1143'}; D.set_par(vCurPar, D.get_db_value('RFR_MRBLOBDATA_PIIM_MR_ACTUAL_BLOB_SIZE') ,'static');
                    vCurPar = {'0x2005', '0x1141'}; D.set_par(vCurPar, D.get_db_value('RFR_MRBLOBDATA_PIIM_MR_COMMENT_STR') ,'static');
                    vCurPar = {'0x2005', '0x1140'}; D.set_par(vCurPar, D.get_db_value('RFR_MRBLOBDATA_PIIM_MR_VERSION_STR') ,'static');
                    vCurPar = {'0x2005', '0x1139'}; D.set_par(vCurPar, D.get_db_value('RFR_MRBLOBDATA_PIIM_MR_TYPE_NAME') ,'static');
                    vCurPar = {'0x2005', '0x1138'}; D.set_par(vCurPar, D.get_db_value('RFR_MRBLOBDATA_PIIM_MR_APPLICATION_NAME') ,'static');
                    vCurPar = {'0x2005', '0x1137'}; D.set_par(vCurPar, D.get_db_value('RFR_MRBLOBDATA_PIIM_MR_BLOB_NAME') ,'static');
                    vCurPar = {'0x2005', '0x1391'}; D.set_par(vCurPar, D.get_db_value('RFR_PATIENT_PIIM_PATIENT_NAME_JOB_IN_PARAMS') ,'static');
                    vCurPar = {'0x0010', '0x4000'}; D.set_par(vCurPar, D.get_db_value('RFR_PATIENT_DICOM_PATIENT_COMMENTS') ,'static');
                    vCurPar = {'0x0010', '0x2160'}; D.set_par(vCurPar, D.get_db_value('RFR_PATIENT_DICOM_ETHNIC_GROUP') ,'static');
                    vCurPar = {'0x0010', '0x0040'}; D.set_par(vCurPar, D.get_db_value('RFR_PATIENT_DICOM_PATIENT_SEX') ,'static');
                    vCurPar = {'0x0010', '0x0030'}; D.set_par(vCurPar, D.get_db_value('RFR_PATIENT_DICOM_PATIENT_BIRTH_DATE', true) ,'static');
                    vCurPar = {'0x0010', '0x0020'}; D.set_par(vCurPar, D.get_db_value('RFR_PATIENT_DICOM_PATIENT_ID', true) ,'static');
                    vCurPar = {'0x0010', '0x0010'}; D.set_par(vCurPar, D.get_db_value('RFR_PATIENT_DICOM_PATIENT_NAME') ,'static');
                    vCurPar = {'0x2050', '0x0020'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_PRESENTATION_LUT_SHAPE') ,'static');
                    vCurPar = {'0x2005', '0x1492'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_SED_VALUE') ,'static');
                    vCurPar = {'0x2005', '0x1480'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_IMAGE_CACHE') ,'static');
                    vCurPar = {'0x2005', '0x144E'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_IS_B1_SERIES') ,'static');
                    vCurPar = {'0x2005', '0x144D'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_IS_B0_SERIES') ,'static');
                    vCurPar = {'0x2005', '0x1448'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_COIL_Q') ,'static');
                    vCurPar = {'0x2005', '0x1447'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_POWER_OPTIMIZATION') ,'static');
                    vCurPar = {'0x2005', '0x1446'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_FWHM_SHIM') ,'static');
                    vCurPar = {'0x2005', '0x1445'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_ATTENUATION_CORRECTION') ,'static');
                    vCurPar = {'0x2005', '0x1444'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_TFE_FACTOR') ,'static');
                    vCurPar = {'0x2005', '0x1443'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_AIM_DDB_DT_LIMIT') ,'static');
                    vCurPar = {'0x2005', '0x1442'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_AIMD_B1_RMS_LIMIT') ,'static');
                    vCurPar = {'0x2005', '0x1441'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_AIMD_WHOLE_BODY_SAR_LIMIT') ,'static');
                    vCurPar = {'0x2005', '0x1440'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_AIMD_HEAD_SAR_LIMIT') ,'static');
                    vCurPar = {'0x2005', '0x143F'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_AIMD_LIMITS_APPLIED') ,'static');
                    vCurPar = {'0x2005', '0x143B'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_IS_COIL_SURVEY') ,'static');
                    vCurPar = {'0x2005', '0x143A'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_DATA_DICTIONARY_CONTENTS_VERSION') ,'static');
                    vCurPar = {'0x2005', '0x1437'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_COLOR_LUT_TYPE') ,'static');
                    vCurPar = {'0x2005', '0x1435'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_SPECTRO_EXAMCARD') ,'static');
                    vCurPar = {'0x2005', '0x1432'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_SAFETY_OVERRIDE_MODE') ,'static');
                    vCurPar = {'0x2005', '0x1431'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_LOCAL_SAR') ,'static');
                    vCurPar = {'0x2005', '0x1414'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_NR_OF_DIFF_B_VALUES') ,'static');
                    vCurPar = {'0x2005', '0x1410'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MF_CONV_TREAT_SPECTOR_MIX_NO') ,'static');
                    vCurPar = {'0x2005', '0x1400'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_VOLUME_VIEW_ENABLED') ,'static');
                    vCurPar = {'0x2005', '0x1399'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MRI_VIEW_BOLD_ENABLED') ,'static');
                    vCurPar = {'0x2005', '0x1398'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_MOBIVIEW_ENABLED') ,'static');
                    vCurPar = {'0x2005', '0x1397'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_ANATOMIC_REG_CODE_VALUE') ,'static');
                    vCurPar = {'0x2005', '0x1396'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_FLOW_IMAGES_PRESENT') ,'static');
                    vCurPar = {'0x2005', '0x1393'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_STATION_NO') ,'static');
                    vCurPar = {'0x2005', '0x1392'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_GEOLINK_ID') ,'static');
                    vCurPar = {'0x2005', '0x1381'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SCANO_GRAM_SURVEY_NUMBER_OF_IMAGES') ,'static');
                    vCurPar = {'0x2005', '0x1370'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_NO_MIXES_SPECTRO') ,'static');
                    vCurPar = {'0x2005', '0x1364'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_VOLUME_SELECTION') ,'static');
                    vCurPar = {'0x2005', '0x1363'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SPECTRUM_PITCH') ,'static');
                    vCurPar = {'0x2005', '0x1362'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SPECTRO_OFFSET') ,'static');
                    vCurPar = {'0x2005', '0x1360'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SPECTRO_VERTICAL_SHIFT') ,'static');
                    vCurPar = {'0x2005', '0x1359'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SPECTRO_TURBO_ECHO_SPACING') ,'static');
                    vCurPar = {'0x2005', '0x1358'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SPECTRO_TITLE_LINE') ,'static');
                    vCurPar = {'0x2005', '0x1357'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SPECTRO_SPECTRAL_BW') ,'static');
                    vCurPar = {'0x2005', '0x1356'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SPECTRO_SI_MODE') ,'static');
                    vCurPar = {'0x2005', '0x1355'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SPECTRO_SICS_INTERVALS') ,'static');
                    vCurPar = {'0x2005', '0x1354'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SPECTRO_SCAN_TYPE') ,'static');
                    vCurPar = {'0x2005', '0x1353'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SPECTRO_PROCESSING_HISTORY') ,'static');
                    vCurPar = {'0x2005', '0x1334'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_PHYSICAL_QUANTITY_FOR_CHEMICAL_SHIFT') ,'static');
                    vCurPar = {'0x2005', '0x1331'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_NUMBER_OF_SPECTRA_ACQUIRED') ,'static');
                    vCurPar = {'0x2005', '0x1329'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SPECTRO_ECHO_TOP_POSITION') ,'static');
                    vCurPar = {'0x2005', '0x1328'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SPECTRO_DATA_ORIGIN') ,'static');
                    vCurPar = {'0x2005', '0x1327'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SPECTRO_COMPLEX_COMPONENT') ,'static');
                    vCurPar = {'0x2005', '0x1326'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SPECTRO_B0_ECHO_TOP_POSITION') ,'static');
                    vCurPar = {'0x2005', '0x1325'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SPECTRO_SI_B0_CORRECTION') ,'static');
                    vCurPar = {'0x2005', '0x1256'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_NR_SC_SOFTWARE_VERSIONS') ,'static');
                    vCurPar = {'0x2005', '0x1249'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_NR_OF_SERIES_OPERATORS_NAME') ,'static');
                    vCurPar = {'0x2005', '0x1245'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_NR_OF_SOFTWARE_VERSION') ,'static');
                    vCurPar = {'0x2005', '0x1134'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_SERIES_TRANSACTION_UID') ,'static');
                    vCurPar = {'0x2005', '0x10C0'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_SERIES_SCAN_SEQUENCE') ,'static');
                    vCurPar = {'0x2005', '0x10A9'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_GEOMETRY_CORRECTION') ,'static');
                    vCurPar = {'0x2005', '0x10A2'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_IS_COCA') ,'static');
                    vCurPar = {'0x2005', '0x109F'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_SPECTRAL_SELECTIVE_EXCITATION_PULSE') ,'static');
                    vCurPar = {'0x2005', '0x1086'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_NUMBER_OF_GEOMETRY') ,'static');
                    vCurPar = {'0x2005', '0x1081'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_STACK_VIEW_AXIS') ,'static');
                    vCurPar = {'0x2005', '0x1070'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_HARD_COPY_PROTOCOL_EV') ,'static');
                    vCurPar = {'0x2005', '0x106F'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_ACQUISITION_TYPE_PRIVATE') ,'static');
                    vCurPar = {'0x2005', '0x103E'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MRRR_INTERVALS_DISTRIBUTION') ,'static');
                    vCurPar = {'0x2005', '0x103D'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_NUMBER_OF_RR_INTERVAL_RANGES') ,'static');
                    vCurPar = {'0x2005', '0x103C'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_TONE') ,'static');
                    vCurPar = {'0x2005', '0x103B'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_TIME_REVERSED_STEADY_STATE') ,'static');
                    vCurPar = {'0x2005', '0x102A'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_PATIENT_REFERENCE_ID') ,'static');
                    vCurPar = {'0x2005', '0x1029'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_PARTIAL_FOURIER_PHASE') ,'static');
                    vCurPar = {'0x2005', '0x1028'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_PARTIAL_FOURIER_FREQUENCY') ,'static');
                    vCurPar = {'0x2005', '0x1027'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_PACKAGE_MODE') ,'static');
                    vCurPar = {'0x2005', '0x1026'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_OVER_SAMPLEING_PHASE') ,'static');
                    vCurPar = {'0x2005', '0x1025'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_NUMBER_OF_VOLUMES') ,'static');
                    vCurPar = {'0x2005', '0x1023'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_NUMBER_OF_SLABS') ,'static');
                    vCurPar = {'0x2005', '0x1022'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_NUMBER_OF_REFERENCES') ,'static');
                    vCurPar = {'0x2005', '0x1021'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_NUMBER_OF_MIXES') ,'static');
                    vCurPar = {'0x2005', '0x1020'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_NR_OF_CHEMICAL_SHIFTS') ,'static');
                    vCurPar = {'0x2005', '0x101F'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_SLICE_SERIES_MPR_PROTOCOL') ,'static');
                    vCurPar = {'0x2005', '0x101E'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_SLICE_SERIES_MIP_PROTOCOL') ,'static');
                    vCurPar = {'0x2005', '0x101D'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_MEASUREMENT_SCAN_RESOLUTION') ,'static');
                    vCurPar = {'0x2005', '0x101C'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_MAGNET_TRANSFER_CONST') ,'static');
                    vCurPar = {'0x2005', '0x101B'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_MAGNETI_PREPARED') ,'static');
                    vCurPar = {'0x2005', '0x101A'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_LABEL_SYNTAX') ,'static');
                    vCurPar = {'0x2005', '0x1019'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_INVERSE_RECONSTRUCTED') ,'static');
                    vCurPar = {'0x2005', '0x1018'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_HARDCOPY_PROTOCOL') ,'static');
                    vCurPar = {'0x2005', '0x1017'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_FOURIER_INTERPOLATION') ,'static');
                    vCurPar = {'0x2005', '0x1016'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_FLOW_COMPENSATION') ,'static');
                    vCurPar = {'0x2005', '0x1015'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_FAT_SATURATION') ,'static');
                    vCurPar = {'0x2005', '0x1014'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_DIFFUSION') ,'static');
                    vCurPar = {'0x2005', '0x1013'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_DEVELOPMENT_MODE') ,'static');
                    vCurPar = {'0x2005', '0x1012'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_CARDIAC_GATING') ,'static');
                    vCurPar = {'0x2001', '0x1060'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_NR_OF_STACKS') ,'static');
                    vCurPar = {'0x2001', '0x105F'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_STACK_SEQUENCE') ,'static');
                    vCurPar = {'0x2001', '0x1025'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_ECHO_TIME_DISPLAY') ,'static');
                    vCurPar = {'0x2001', '0x1024'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_IS_INTERACTIVE') ,'static');
                    vCurPar = {'0x2001', '0x1023'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_FLIP_ANGLE') ,'static');
                    vCurPar = {'0x2001', '0x1022'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_WATER_FAT_SHIFT') ,'static');
                    vCurPar = {'0x2001', '0x1021'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_SEL_PART_INVERSION_RECOVERY') ,'static');
                    vCurPar = {'0x2001', '0x1020'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_SCANNING_TECHNIQUE_DESC') ,'static');
                    vCurPar = {'0x2001', '0x101F'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_RESPIRATION_SYNC') ,'static');
                    vCurPar = {'0x2001', '0x101E'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_REFORMAT_ACCURACY') ,'static');
                    vCurPar = {'0x2001', '0x101C'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_PREPULSE_TYPE') ,'static');
                    vCurPar = {'0x2001', '0x101B'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_PREPULSE_DELAY') ,'static');
                    vCurPar = {'0x2001', '0x1019'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_PARTIAL_MATRIX_SCANNED') ,'static');
                    vCurPar = {'0x2001', '0x1018'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_NR_OF_SLICES') ,'static');
                    vCurPar = {'0x2001', '0x1017'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_NR_OF_PHASES') ,'static');
                    vCurPar = {'0x2001', '0x1016'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_NR_OF_PHASE_CONTRAST_DIRCTNS') ,'static');
                    vCurPar = {'0x2001', '0x1015'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_NR_OF_LOCATIONS') ,'static');
                    vCurPar = {'0x2001', '0x1014'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_NR_OF_ECHOES') ,'static');
                    vCurPar = {'0x2001', '0x1013'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_EPI_FACTOR') ,'static');
                    vCurPar = {'0x2001', '0x1012'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_DYNAMIC_SERIES') ,'static');
                    vCurPar = {'0x0028', '0x0010'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_ROWS') ,'static');
                    vCurPar = {'0x0020', '0x9256'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_RESPIRATORY_TRIGGER_DELAY_TRESHOLD') ,'static');
                    vCurPar = {'0x0020', '0x9255'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_RESPIRATORY_TRIGGER_DELAY_TIME') ,'static');
                    vCurPar = {'0x0020', '0x9254'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_RESPIRATORY_INTERVAL_TIME') ,'static');
                    vCurPar = {'0x0020', '0x9072'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_FRAME_LATERALITY') ,'static');
                    vCurPar = {'0x0020', '0x1040'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_POSITION_REFERENCE_INDICATOR') ,'static');
                    vCurPar = {'0x0020', '0x0060'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_LATERALITY') ,'static');
                    vCurPar = {'0x0020', '0x0052'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_FRAME_OF_REFERENCE_UID') ,'static');
                    vCurPar = {'0x0018', '0x9241'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_GRADIENT_ECHO_TRAIN_LENGTH') ,'static');
                    vCurPar = {'0x0018', '0x9240'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_RF_ECHO_TRAIN_LENGTH') ,'static');
                    vCurPar = {'0x0018', '0x9231'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_MR_ACQUISITION_PHASE_ENCODING_STEPS_IN_PLANE') ,'static');
                    vCurPar = {'0x0018', '0x9200'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_MR_SPECTROSCOPY_ACQUISITION_TYPE') ,'static');
                    vCurPar = {'0x0018', '0x9199'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_WATER_REFERENCED_PHASE_CORRECTION') ,'static');
                    vCurPar = {'0x0018', '0x9183'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_FLOW_COMPENSATION_DIRECTION') ,'static');
                    vCurPar = {'0x0018', '0x9182'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_GRADIENT_OUTPUT') ,'static');
                    vCurPar = {'0x0018', '0x9181'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SPECIFIC_ABSORPTION_RATE_VALUE') ,'static');
                    vCurPar = {'0x0018', '0x9180'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_GRADIENT_OUTPUT_TYPE') ,'static');
                    vCurPar = {'0x0018', '0x9179'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SPECIFIC_ABSORPTION_RATE_DEFINITION') ,'static');
                    vCurPar = {'0x0018', '0x9174'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_APPLICABLE_SAFETY_STANDARD_AGENCY') ,'static');
                    vCurPar = {'0x0018', '0x9172'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_BULK_MOTION_COMPENSATION_TECHNIQUE') ,'static');
                    vCurPar = {'0x0018', '0x9171'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_RESPIRATORY_SIGNAL_SOURCE') ,'static');
                    vCurPar = {'0x0018', '0x9170'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_RESPIRATORY_MOTION_COMPENSATION_TECHNIQUE') ,'static');
                    vCurPar = {'0x0018', '0x9060'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_DE_COUPLED_NUCLEUS') ,'static');
                    vCurPar = {'0x0018', '0x9059'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_DE_COUPLING') ,'static');
                    vCurPar = {'0x0018', '0x9058'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_MR_ACQUISITION_FREQUENCY_ENCODING_STEPS') ,'static');
                    vCurPar = {'0x0018', '0x9053'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_CHEMICAL_SHIFT_REFERENCE') ,'static');
                    vCurPar = {'0x0018', '0x9051'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_TRANSMIT_COIL_TYPE') ,'static');
                    vCurPar = {'0x0018', '0x9050'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_TRANSMIT_COIL_MANUFACTURER_NAME') ,'static');
                    vCurPar = {'0x0018', '0x9037'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_CARDIAC_SYNCHRONIZATION_TECHNIQUE') ,'static');
                    vCurPar = {'0x0018', '0x9036'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_PARTIAL_FOURIER_DIRECTION') ,'static');
                    vCurPar = {'0x0018', '0x9035'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_TAG_THICKNESS') ,'static');
                    vCurPar = {'0x0018', '0x9034'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_RECTILINEAR_PHASE_ENCODE_REORDERING') ,'static');
                    vCurPar = {'0x0018', '0x9033'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SEGMENTED_K_SPACE_TRAVERSAL') ,'static');
                    vCurPar = {'0x0018', '0x9032'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_GEOMETRY_OF_K_SPACE_TRAVERSAL') ,'static');
                    vCurPar = {'0x0018', '0x9030'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_TAG_SPACING_FIRST_DIMENSION') ,'static');
                    vCurPar = {'0x0018', '0x9029'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_OVERSAMPLING_PHASE') ,'static');
                    vCurPar = {'0x0018', '0x9028'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_TAGGING') ,'static');
                    vCurPar = {'0x0018', '0x9026'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SPECTRALLY_SELECTED_EXCITATION') ,'static');
                    vCurPar = {'0x0018', '0x9025'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SPECTRALLY_SELECTED_SUPPRESSION') ,'static');
                    vCurPar = {'0x0018', '0x9024'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SATURATION_RECOVERY') ,'static');
                    vCurPar = {'0x0018', '0x9022'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_BLOOD_SIGNAL_NULLING') ,'static');
                    vCurPar = {'0x0018', '0x9021'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_T2_PREPARATION') ,'static');
                    vCurPar = {'0x0018', '0x9020'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_MAGNETIZATION_TRANSFER') ,'static');
                    vCurPar = {'0x0018', '0x9019'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_TAG_ANGLE_FIRST_AXIS') ,'static');
                    vCurPar = {'0x0018', '0x1016'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SECONDARY_CAPTURE_DEVICE_MANUFACTURER') ,'static');
                    vCurPar = {'0x0018', '0x1010'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SECONDARY_CAPTURE_DEVICE_ID') ,'static');
                    vCurPar = {'0x0018', '0x1000'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_DEVICE_SERIAL_NUMBER', true) ,'static');
                    vCurPar = {'0x0018', '0x0095'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_PIXEL_BANDWIDTH') ,'static');
                    vCurPar = {'0x0018', '0x0094'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_PERCENT_PHASE_FIELD_OF_VIEW') ,'static');
                    vCurPar = {'0x0018', '0x0093'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_PERCENT_SAMPLING') ,'static');
                    vCurPar = {'0x0018', '0x0091'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_ECHO_TRAIN_LENGTH') ,'static');
                    vCurPar = {'0x0018', '0x0089'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_NUMBER_OF_PHASE_ENCODING_STEPS') ,'static');
                    vCurPar = {'0x0018', '0x0082'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_INVERSION_TIME') ,'static');
                    vCurPar = {'0x0018', '0x0015'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_BODY_PART_EXAMINED') ,'static');
                    vCurPar = {'0x0008', '0x9207'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_VOLUME_BASED_CALCULATION_TECHNIQUE') ,'static');
                    vCurPar = {'0x0008', '0x9206'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_VOLUMETRIC_PROPERTIES') ,'static');
                    vCurPar = {'0x0008', '0x9205'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_PIXEL_PRESENTATION') ,'static');
                    vCurPar = {'0x0008', '0x9123'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_CREATOR_VERSION_UID') ,'static');
                    vCurPar = {'0x0008', '0x1090'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_MANUFACTURERS_MODEL_NAME') ,'static');
                    vCurPar = {'0x0008', '0x1070'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_OPERATORS_NAME') ,'static');
                    vCurPar = {'0x0008', '0x1050'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_PERFORMING_PHYSICIANS_NAME') ,'static');
                    vCurPar = {'0x0008', '0x1040'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_INSTITUTIONAL_DEPARTMENT_NAME') ,'static');
                    vCurPar = {'0x0008', '0x103E'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SERIES_DESCRIPTION') ,'static');
                    vCurPar = {'0x0008', '0x1010'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_STATION_NAME') ,'static');
                    vCurPar = {'0x0008', '0x0104'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_CODE_MEANING') ,'static');
                    vCurPar = {'0x0008', '0x0102'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_CODING_SCHEME_DESIGNATOR') ,'static');
                    vCurPar = {'0x0008', '0x0100'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_CODE_VALUE') ,'static');
                    vCurPar = {'0x0008', '0x0081'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_INSTITUTION_ADDRESS') ,'static');
                    vCurPar = {'0x2005', '0x1406'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIESBLOBSET_PIIM_MR_NR_OF_SPECIFIC_CHARACTER_SET') ,'static');
                    vCurPar = {'0x2005', '0x1132'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIESBLOBSET_PIIM_MR_BLOB_DATA_OBJECT_ARRAY_SEQUENCE') ,'static');
                    vCurPar = {'0x0008', '0x0014'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIESBLOBSET_DICOM_INSTANCE_CREATOR_UID') ,'static');
                    vCurPar = {'0x0008', '0x0005'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIESBLOBSET_DICOM_SPECIFIC_CHARACTER_SET') ,'static');
                    vCurPar = {'0x2005', '0x143E'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_MR_STACK_POSTERIOR_COIL_POS') ,'static');
                    vCurPar = {'0x2005', '0x143D'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_MR_STACK_TABLE_POS_LAT') ,'static');
                    vCurPar = {'0x2005', '0x143C'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_MR_STACK_TABLE_POS_LONG') ,'static');
                    vCurPar = {'0x2005', '0x1081'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_MR_STACK_VIEW_AXIS') ,'static');
                    vCurPar = {'0x2005', '0x107E'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_MR_STACK_SLICE_DISTANCE') ,'static');
                    vCurPar = {'0x2005', '0x107B'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_MR_STACK_PREPARATION_DIRECTION') ,'static');
                    vCurPar = {'0x2005', '0x107A'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_MR_STACK_OFFCENTRE_RL') ,'static');
                    vCurPar = {'0x2005', '0x1079'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_MR_STACK_OFFCENTRE_FH') ,'static');
                    vCurPar = {'0x2005', '0x1078'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_MR_STACK_OFFCENTRE_AP') ,'static');
                    vCurPar = {'0x2005', '0x1076'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_MR_STACK_FOV_RL') ,'static');
                    vCurPar = {'0x2005', '0x1075'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_MR_STACK_FOV_FH') ,'static');
                    vCurPar = {'0x2005', '0x1074'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_MR_STACK_FOV_AP') ,'static');
                    vCurPar = {'0x2005', '0x1073'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_MR_STACK_ANGULATION_RL') ,'static');
                    vCurPar = {'0x2005', '0x1072'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_MR_STACK_ANGULATION_FH') ,'static');
                    vCurPar = {'0x2005', '0x1071'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_MR_STACK_ANGULATION_AP') ,'static');
                    vCurPar = {'0x2001', '0x1036'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_STACK_TYPE') ,'static');
                    vCurPar = {'0x2001', '0x1035'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_STACK_SLICE_NUMBER') ,'static');
                    vCurPar = {'0x2001', '0x1033'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_STACK_RADIAL_AXIS') ,'static');
                    vCurPar = {'0x2001', '0x1032'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_STACK_RADIAL_ANGLE') ,'static');
                    vCurPar = {'0x2001', '0x102D'}; D.set_par(vCurPar, D.get_db_value('RFR_STACK_PIIM_STACK_NUMBER_OF_SLICES') ,'static');
                    vCurPar = {'0x2005', '0x1401'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_PIIM_MR_NUMBER_OF_STUDY_REFERENCE') ,'static');
                    vCurPar = {'0x2005', '0x1382'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_PIIM_MR_NUMBER_OF_PROCEDURE_CODES') ,'static');
                    vCurPar = {'0x2005', '0x1253'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_PIIM_MR_NR_OF_STUDY_PATIENT_MEDICAL_ALERTS') ,'static');
                    vCurPar = {'0x2005', '0x1252'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_PIIM_MR_NR_OF_STUDY_PATIENT_CONTRAST_ALLERGIES') ,'static');
                    vCurPar = {'0x2005', '0x1251'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_PIIM_MR_NR_OF_STUDY_ADMITTING_DIAGNOSTIC_DESCR') ,'static');
                    vCurPar = {'0x2005', '0x1060'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_PIIM_MR_STUDY_SEQUENCE_NUMBER') ,'static');
                    vCurPar = {'0x2005', '0x105F'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_PIIM_MR_STUDY_ORIGIN') ,'static');
                    vCurPar = {'0x0040', '0x2400'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_IMAGING_SERVICE_REQUEST_COMMENTS') ,'static');
                    vCurPar = {'0x0040', '0x2010'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_ORDER_CALLBACK_PHONE_NUMBER') ,'static');
                    vCurPar = {'0x0040', '0x2009'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_ORDER_ENTERERS_LOCATION') ,'static');
%                     vCurPar = {'0x0040', '0x2005'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_ISSUE_TIME_OF_IMAGING_SERVICE_REQUEST') ,'static');
%                     vCurPar = {'0x0040', '0x2004'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_ISSUE_DATE_OF_IMAGING_SERVICE_REQUEST') ,'static');
                    vCurPar = {'0x0040', '0x2001'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_REASON_FOR_THE_IMAGING_SERVICE_REQUEST') ,'static');
                    vCurPar = {'0x0040', '0x1400'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_REQUESTED_PROCEDURE_COMMENTS') ,'static');
                    vCurPar = {'0x0040', '0x1005'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_REQUESTED_PROCEDURE_LOCATION') ,'static');
                    vCurPar = {'0x0040', '0x1004'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_PATIENT_TRANSPORT_ARRANGEMENTS') ,'static');
                    vCurPar = {'0x0040', '0x1003'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_REQUESTED_PROCEDURE_PRIORITY') ,'static');
                    vCurPar = {'0x0040', '0x1002'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_REASON_FOR_THE_REQUESTED_PROCEDURE') ,'static');
                    vCurPar = {'0x0040', '0x1001'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_REQUESTED_PROCEDURE_ID', true) ,'static'); % LMG debug: force char output                    
                    vCurPar = {'0x0038', '0x0500'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_PATIENT_STATE') ,'static');
                    vCurPar = {'0x0038', '0x0050'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_SPECIAL_NEEDS') ,'static');
                    vCurPar = {'0x0032', '0x4000'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_MR_STUDY_COMMENTS') ,'static');
                    vCurPar = {'0x0008', '0x0030'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_STUDY_TIME', true) ,'static');
                    vCurPar = {'0x0008', '0x0020'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_STUDY_DATE', true) ,'static');
                    
                    
                    % manually added with the help of "find_rfr_pars_for_dicom_tags.m" and "choose_dicom_options.m"
                    vCurPar = 'SeriesTime'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SERIES_TIME', true) ,'static');
                    vCurPar = 'AccessionNumber'; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_ACCESSION_NUMBER', true) ,'static');
                    vCurPar = 'Modality'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_MODALITY') ,'static');
                    vCurPar = 'Manufacturer'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_MANUFACTURER') ,'static');
                    vCurPar = 'InstitutionName'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_INSTITUTION_NAME') ,'static');
%                     vCurPar = 'CodeValue'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_ANATOMIC_REG_CODE_VALUE') ,'static');
                    vCurPar = 'StudyDescription'; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_STUDY_DESCRIPTION') ,'static');
                    vCurPar = 'OperatorsName'; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_DICOM_OPERATORS_NAME') ,'static');
                    vCurPar = 'AdmittingDiagnosesDescription'; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_ADMITTING_DIAGNOSES_DESCRIPTION') ,'static');
                    vCurPar = 'OtherPatientIDs'; D.set_par(vCurPar, D.get_db_value('RFR_PATIENT_DICOM_OTHER_PATIENT_IDS') ,'static');
                    vCurPar = 'PatientWeight'; D.set_par(vCurPar, D.get_db_value('EX_ACQ_patient_weight') ,'static');
                    vCurPar = 'PregnancyStatus'; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_PREGNANCY_STATUS') ,'static');
                    vCurPar = 'BodyPartExamined'; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_DICOM_BODY_PART_EXAMINED') ,'static');
                    vCurPar = 'SliceThickness'; D.set_par(vCurPar, D.get_db_value('RC_slice_thickness') ,'static');
                    vCurPar = 'RepetitionTime'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_REPETITION_TIME') ,'static');
                    vCurPar = 'NumberOfAverages'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_NUMBER_OF_AVERAGES') ,'static');
                    vCurPar = 'ImagingFrequency'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_IMAGING_FREQUENCY') ,'static');
                    vCurPar = 'ImagedNucleus'; D.set_par(vCurPar, D.get_db_value('VW_imaged_nucleus') ,'static');
                    vCurPar = 'MagneticFieldStrength'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_MAGNETIC_FIELD_STRENGTH') ,'static');
                    vCurPar = 'EchoTrainLength'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_ECHO_TRAIN_LENGTH') ,'static');
                    vCurPar = 'PercentSampling'; D.set_par(vCurPar, D.get_db_value('VW_percent_sampling') ,'static');
                    vCurPar = 'SoftwareVersions'; D.set_par(vCurPar, D.get_db_value({'RFR_SERIES_DICOM_SOFTWARE_VERSIONS0', 'RFR_SERIES_DICOM_SOFTWARE_VERSIONS1'}) ,'static');
                    vCurPar = 'SoftwareVersion'; D.set_par(vCurPar, D.get_db_value({'RFR_SERIES_DICOM_SOFTWARE_VERSIONS0', 'RFR_SERIES_DICOM_SOFTWARE_VERSIONS1'}) ,'static'); % LMG edit/bugfix
                    vCurPar = 'ProtocolName'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_PROTOCOL_NAME') ,'static');
                    vCurPar = 'HeartRate'; D.set_par(vCurPar, D.get_db_value('VW_cardiac_freq') ,'static');
                    vCurPar = 'FlipAngle'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_FLIP_ANGLE') ,'static');
                    vCurPar = 'SAR'; D.set_par(vCurPar, D.get_db_value('AC_act_SAR') ,'static');
                    vCurPar = 'dBdt'; D.set_par(vCurPar, D.get_db_value('AC_dbdt_level') ,'static');
                    vCurPar = 'StudyInstanceUID'; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_STUDY_INSTANCE_UID') ,'static');
                    
                    % Calculate the Study ID. The Study ID is a short string and
                    % cannot be longer than 16 chars. Therefore try to extract the
                    % date from the StudyInstanceUID
                    study_instance_uid = D.get_db_value('RFR_STUDY_DICOM_STUDY_INSTANCE_UID');
                    if length(study_instance_uid) >= 16;
                        study_id = study_instance_uid(end-15:end);
                    else
                        study_id = study_instance_uid;
                    end
                    dotind = strfind(study_instance_uid,'.');
                    if( ~isempty(dotind) && dotind(end) < length(study_instance_uid))
                        dotind = dotind(end);
                        date = study_instance_uid(dotind+1:end);
                        try
                            DateNumber = datenum(date, 'yyyymmddHHMMSSFFF');
                            study_id = round(DateNumber*10000);
                        catch
                        end                    
                    end                                           
                    vCurPar = 'StudyID'; D.set_par(vCurPar, num2str(study_id) ,'static');
                    
                    vCurPar = 'AcquisitionNumber'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_ACQUISITION_NUMBER') ,'static');
                    vCurPar = 'Columns'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SIGNAL_DOMAIN_COLUMNS') ,'static');
                    vCurPar = 'PerformedStationAETitle'; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_PERFORMED_STATION_AE_TITLE') ,'static');
                    vCurPar = 'PerformedProcedureStepStartDate'; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_PERFORMED_PROCEDURE_STEP_START_DATE', true) ,'static');
                    vCurPar = 'PerformedProcedureStepStartTime'; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_PERFORMED_PROCEDURE_STEP_START_TIME', true) ,'static');
                    vCurPar = 'PerformedProcedureStepEndDate'; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_PERFORMED_PROCEDURE_STEP_END_DATE', true) ,'static');
                    vCurPar = 'PerformedProcedureStepEndTime'; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_PERFORMED_PROCEDURE_STEP_END_TIME', true) ,'static');
                    vCurPar = 'PerformedProcedureStepID'; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_PERFORMED_PROCEDURE_STEP_ID', true) ,'static');
                    vCurPar = 'PerformedProcedureStepDescription'; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_PERFORMED_PROCEDURE_STEP_DESCRIPTION') ,'static');
%                     vCurPar = 'CodeValue'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_ANATOMIC_REG_CODE_VALUE') ,'static');
                    vCurPar = {'0x2001', '0x100c'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_ARRHYTHMIA_REJECTION') ,'static');
                    vCurPar = {'0x2001', '0x100e'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_CARDIAC_CYCLED') ,'static');
                    vCurPar = {'0x2001', '0x100f'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_CARDIAC_GATE_WIDTH') ,'static');
                    vCurPar = {'0x2001', '0x1010'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_CARDIAC_SYNC') ,'static');
                    vCurPar = {'0x2001', '0x1011'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_DIFFUSION_ECHO_TIME') ,'static');
                    vCurPar = {'0x2001', '0x1032'}; D.set_par(vCurPar, D.get_db_value('EX_GEO_cur_stack_radial_angle') ,'static');
                    vCurPar = {'0x2001', '0x1033'}; D.set_par(vCurPar, D.get_db_value('EX_GEO_cur_stack_radial_axis') ,'static');
                    vCurPar = {'0x2001', '0x1036'}; D.set_par(vCurPar, upper(D.get_db_value('EX_GEO_cur_stack_type')) ,'static');
                    vCurPar = {'0x2001', '0x1063'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_PIIM_EXAMINATION_SOURCE') ,'static');
                    if strcmp(D.mRelease, 'REL53')
                        vCurPar = {'0x2005', '0x101e'}; D.set_par(vCurPar, D.get_db_value('VAL01_PV_mip_protocol') ,'static');
                        vCurPar = {'0x2005', '0x101f'}; D.set_par(vCurPar, D.get_db_value('VAL01_PV_mip_protocol') ,'static');
                    else
                        vCurPar = {'0x2005', '0x101e'}; D.set_par(vCurPar, D.get_db_value('UGN1_PV_mip_protocol') ,'static');
                        vCurPar = {'0x2005', '0x101f'}; D.set_par(vCurPar, D.get_db_value('UGN1_PV_mpr_protocol') ,'static');
                    end
%                     vCurPar = 'CodeValue'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_ANATOMIC_REG_CODE_VALUE') ,'static');
                    vCurPar = 'Spoiling'; D.set_par(vCurPar, D.get_db_value('CSC_rear_spoiling') ,'static');
                    vCurPar = 'Tagging'; D.set_par(vCurPar, D.get_db_value('IF_str_ctagging') ,'static');
                    vCurPar = 'ChemicalShiftReference'; D.set_par(vCurPar, D.get_db_value({'RFR_SERIES_DICOM_CHEMICAL_SHIFT_REFERENCE0', 'RFR_SERIES_DICOM_CHEMICAL_SHIFT_REFERENCE1'}) ,'static');
                    vCurPar = 'ParallelReductionFactorInPlane'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_PARALLEL_REDUCTION_FACTOR_IN_PLANE') ,'static');
                    vCurPar = 'ParallelAcquisition'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_PARALLEL_ACQUISITION') ,'static');
                    vCurPar = 'PartialFourier'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_PARTIAL_FOURIER') ,'static');
                    vCurPar = 'VelocityEncodingDirection'; D.set_par(vCurPar, D.get_db_value({'RFR_SERIES_DICOM_VELOCITY_ENCODING_DIRECTION0', 'RFR_SERIES_DICOM_VELOCITY_ENCODING_DIRECTION1', 'RFR_SERIES_DICOM_VELOCITY_ENCODING_DIRECTION2'}) ,'static');
                    vCurPar = 'VelocityEncodingMinimumValue'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_VELOCITY_ENCODING_MINIMUM_VALUE') ,'static');
                    vCurPar = 'NumberOfKSpaceTrajectories'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_NUMBER_OF_K_SPACE_TRAJECTORIES') ,'static');
                    vCurPar = 'FrequencyCorrection'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_FREQUENCY_CORRECTION') ,'static');
                    vCurPar = 'ParallelReductionFactorOutOfPlane'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_PARALLEL_REDUCTION_FACTOR_OUT_OF_PLANE') ,'static');
                    vCurPar = 'ParallelReductionFactorSecondInPlane'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_PARALLEL_REDUCTION_FACTOR_SECOND_IN_PLANE') ,'static');
                    vCurPar = 'CardiacBeatRejectionTechnique'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_CARDIAC_BEAT_REJECTION_TECHNIQUE') ,'static');
                    vCurPar = 'LUTExplanation'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_LUT_EXPLANATION') ,'static');
                    vCurPar = 'ConversionType'; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_DICOM_CONVERSION_TYPE') ,'static');
                    vCurPar = 'InstitutionAddress'; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_DICOM_INSTITUTION_ADDRESS') ,'static');
                    vCurPar = 'PatientWeight'; D.set_par(vCurPar, D.get_db_value('EX_ACQ_patient_weight') ,'static');
                    vCurPar = 'MedicalAlerts'; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_MEDICAL_ALERTS') ,'static');
                    vCurPar = 'Allergies'; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_CONTRAST_ALLERGIES') ,'static');
                    vCurPar = 'Occupation'; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_OCCUPATION') ,'static');
                    vCurPar = 'SliceThickness'; D.set_par(vCurPar, D.get_db_value('RC_slice_thickness') ,'static');
                    vCurPar = 'PercentSampling'; D.set_par(vCurPar, D.get_db_value('VW_percent_sampling') ,'static');
                    vCurPar = 'SecondaryCaptureDeviceSoftwareVersions'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SECONDARY_CAPTURE_DEVICE_SOFTWARE_VERSIONS') ,'static');
                    vCurPar = 'VideoImageFormatAcquired'; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_DICOM_VIDEO_IMAGE_FORMAT_ACQUIRED') ,'static');
                    vCurPar = 'DigitalImageFormatAcquired'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_DIGITAL_IMAGE_FORMAT_ACQUIRED') ,'static');
                    vCurPar = 'HeartRate'; D.set_par(vCurPar, D.get_db_value('VW_cardiac_freq') ,'static');
                    vCurPar = 'SAR'; D.set_par(vCurPar, D.get_db_value('AC_act_SAR') ,'static');
                    vCurPar = 'dBdt'; D.set_par(vCurPar, D.get_db_value('AC_dbdt_level') ,'static');
                    vCurPar = 'RequestingPhysician'; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_REQUESTING_PHYSICIAN') ,'static');
                    vCurPar = 'RequestingService'; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_REQUESTING_SERVICE') ,'static');
                    vCurPar = 'RequestedProcedureDescription'; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_REQUESTED_PROCEDURE_DESCRIPTION') ,'static');
                    vCurPar = 'RequestedContrastAgent'; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_REQUESTED_CONTRAST_AGENT') ,'static');
                    vCurPar = 'PerformedStationName'; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_PERFORMED_STATION_NAME') ,'static');
                    vCurPar = 'PerformedLocation'; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_PERFORMED_LOCATION') ,'static');
                    vCurPar = 'PerformedProcedureStepStatus'; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_PERFORMED_PROCEDURE_STEP_STATUS') ,'static');
                    vCurPar = 'PerformedProcedureTypeDescription'; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_DICOM_PERFORMED_PROCEDURE_TYPE_DESCRIPTION') ,'static');
                    vCurPar = 'AcquisitionDateTime'; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_DICOM_SERIES_DATE', true) ,'static');
                    vCurPar = 'DecouplingMethod'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_DE_COUPLING_METHOD') ,'static');
                    vCurPar = 'SignalDomainColumns'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SIGNAL_DOMAIN_COLUMNS') ,'static');
                    vCurPar = 'DataRepresentation'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_DATA_REPRESENTATION') ,'static');
                    
                    % manually added
                    vCurPar = 'AcquisitionDate'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SERIES_DATE', true) ,'static');
                    vCurPar = 'ContentDate'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SERIES_DATE', true) ,'static');
                    vCurPar = 'PerformingPhysicianName'; D.set_par(vCurPar, D.get_db_value('RFR_ECSERIES_DICOM_PERFORMING_PHYSICIANS_NAME') ,'static');
                    vCurPar = 'PatientWeight'; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_DICOM_PATIENTS_WEIGHT') ,'static');
                    vCurPar = 'ScanningSequence'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_SERIES_SCAN_SEQUENCE') ,'static');
                    vCurPar = 'NumberOfPhaseEncodingSteps'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_MR_ACQUISITION_PHASE_ENCODING_STEPS_IN_PLANE') ,'static');
                    vCurPar = 'SecondaryCaptureDeviceManufacturerModelName'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_DICOM_SECONDARY_CAPTURE_DEVICE_MANUFACTURERS_MODEL_NAME') ,'static');
                    vCurPar = 'NumberOfTemporalPositions'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_NR_OF_DYNAMIC_SCANS') ,'static');
                    vCurPar = 'LowRRValue'; D.set_par(vCurPar, D.get_db_value('VW_min_rr_interval') ,'static');
                    vCurPar = 'HighRRValue'; D.set_par(vCurPar, D.get_db_value('VW_max_rr_interval') ,'static');
                    vCurPar = 'IntervalsRejected'; D.set_par(vCurPar, D.get_db_value('VW_num_of_rejections') ,'static');
                    intervals = D.get_db_value('VW_rr_intervals');
                    if ~isempty(intervals)
                        intervals = intervals(1);
                    end
                    vCurPar = 'IntervalsAcquired'; D.set_par(vCurPar, intervals ,'static');
                    vCurPar = {'0x2001', '0x1019'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_PARTIAL_MATRIX_SCANNED') ,'static');
                    vCurPar = {'0x2001', '0x101b'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_PREPULSE_DELAY') ,'static');
                    vCurPar = {'0x2001', '0x101c'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_PREPULSE_TYPE') ,'static');
                    vCurPar = 'PulseSequenceName'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_SCANNING_TECHNIQUE_DESC') ,'static');
                    vCurPar = 'InversionTimes'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_INVERSION_TIME') ,'static');
                    vCurPar = 'InversionTimes'; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_NR_OF_CHEMICAL_SHIFTS') ,'static');
                    vCurPar = {'0x2001', '0x1081'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_NR_OF_DYNAMIC_SCANS') ,'static');
                    vCurPar = {'0x2005', '0x1428'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_NR_OF_LABEL_TYPES') ,'static');
                    vCurPar = {'0x2001', '0x1086'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_NR_OF_PHASE_ENCODING_STEPS') ,'static');
                    
                    
                    % DB parameters for which the value is directly set here
                    if( ~strcmpi(D.get_db_value('RFR_SERIES_PIIM_MR_SERIES_SEL_PART_INVERSION_RECOVERY'), 'N') ) inversion_recovery = 'YES'; else inversion_recovery = 'NO'; end
                    vCurPar = 'InversionRecovery'; D.set_par(vCurPar, inversion_recovery ,'static');
                    if( ~strcmpi(D.get_db_value('RFR_SERIES_PIIM_MR_SPATIAL_PRESATURATION'), 'N') ) presat = 'SLAB'; else presat = 'NONE'; end
                    vCurPar = 'SpatialPresaturation'; D.set_par(vCurPar, presat ,'static');
                else
                    % ---------------------------------------------
                    % Release 3 and higher
                    % ---------------------------------------------
                    vCurPar = {'0x2005', '0x1442'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_AIMD_B1_RMS_limit') ,'static');
                    vCurPar = {'0x2005', '0x1440'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_AIMD_head_SAR_limit') ,'static');
                    vCurPar = {'0x2005', '0x143F'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_AIMD_limits_applied') ,'static');
                    vCurPar = {'0x2005', '0x1441'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_AIMD_wholebody_SAR_lim') ,'static');
                    vCurPar = {'0x2005', '0x1443'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_AIMD_dbDt_limit') ,'static');
                    vCurPar = {'0x0008', '0x0050'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_accession_no') ,'static');
                    vCurPar = {'0x1001', '0x1009'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_series_archived_status') ,'static');
                    vCurPar = {'0x2005', '0x1445'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_attenuation_correction') ,'static');
                    vCurPar = {'0x0018', '0x9022'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_blood_signal_nulling') ,'static');
                    vCurPar = {'0x2005', '0x1448'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_coil_q') ,'static');
                    vCurPar = {'0x2005', '0x1437'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_color_lut_type') ,'static');
                    vCurPar = {'0x2005', '0x1456'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_contrast_no_injections') ,'static');
                    vCurPar = {'0x0018', '0x9094'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_coverage_of_k_space') ,'static');
                    vCurPar = {'0x0018', '0x9008'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_echo_pulse_sequence') ,'static');
                    vCurPar = {'0x2001', '0x10c8'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_card_name') ,'static');
                    vCurPar = {'0x2005', '0x142B'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_export_status') ,'static');
                    vCurPar = {'0x2005', '0x142D'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_mediawrite_status') ,'static');
                    vCurPar = {'0x2005', '0x142A'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_print_status') ,'static');
                    vCurPar = {'0x2005', '0x142C'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_storcommit_status') ,'static');
                    vCurPar = {'0x0040', '0x0280'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_comments') ,'static');
                    vCurPar = {'0x0040', '0x0254'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_desc') ,'static');
                    vCurPar = {'0x0040', '0x0250'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_end_date') ,'static');
                    vCurPar = {'0x0040', '0x0251'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_end_time') ,'static');
                    vCurPar = {'0x0040', '0x0253'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_examination_id') ,'static');
                    vCurPar = {'0x0040', '0x0243'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_location') ,'static');
                    vCurPar = {'0x2001', '0x1063'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_source') ,'static');
                    vCurPar = {'0x0040', '0x0244'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_start_date', true) ,'static');
                    vCurPar = {'0x0040', '0x0245'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_start_time', true) ,'static');
                    vCurPar = {'0x0040', '0x0241'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_station_ae_title') ,'static');
                    vCurPar = {'0x0040', '0x0242'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_station_name') ,'static');
                    vCurPar = {'0x0040', '0x0252'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_status') ,'static');
                    vCurPar = {'0x0040', '0x0255'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_exam_type_desc') ,'static');
                    vCurPar = {'0x1001', '0x1008'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_export_status') ,'static');
                    vCurPar = {'0x2005', '0x1446'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_fwhm_shim') ,'static');
                    vCurPar = {'0x0020', '0x9072'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_frame_laterality') ,'static');
                    vCurPar = {'0x0020', '0x1040'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_position_ref_indicator') ,'static');
                    vCurPar = {'0x0020', '0x0052'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_frame_of_ref_uid') ,'static');
                    vCurPar = {'0x0018', '0x9182'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_gradient_output') ,'static');
                    vCurPar = {'0x2005', '0x144D'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_is_b0_series') ,'static');
                    vCurPar = {'0x2005', '0x144E'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_is_b1_series') ,'static');
                    vCurPar = {'0x2005', '0x1420'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_lut1_begin_color') ,'static');
                    vCurPar = {'0x2005', '0x1421'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_lut1_end_color') ,'static');
                    vCurPar = {'0x2005', '0x141E'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_lut1_offset') ,'static');
                    vCurPar = {'0x2005', '0x141F'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_lut1_range') ,'static');
                    vCurPar = {'0x2005', '0x1424'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_lut2_begin_color') ,'static');
                    vCurPar = {'0x2005', '0x1425'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_lut2_end_color') ,'static');
                    vCurPar = {'0x2005', '0x1422'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_lut2_offset') ,'static');
                    vCurPar = {'0x2005', '0x1423'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_lut2_range') ,'static');
                    vCurPar = {'0x1001', '0x100B'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_media_write_status') ,'static');
                    vCurPar = {'0x0018', '0x9058'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_mr_acq_freq_enc_steps') ,'static');
                    vCurPar = {'0x0018', '0x9231'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_mr_acq_ph_enc_stp_inp') ,'static');
                    vCurPar = {'0x2005', '0x1012'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_cardiac_gating') ,'static');
                    vCurPar = {'0x2005', '0x1015'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_fat_saturation') ,'static');
                    vCurPar = {'0x2005', '0x1016'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_flow_compensation') ,'static');
                    vCurPar = {'0x2005', '0x1396'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_flow_images_present') ,'static');
                    vCurPar = {'0x2005', '0x1017'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_fourier_interpolation') ,'static');
                    vCurPar = {'0x2005', '0x1392'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_geolink_id') ,'static');
                    vCurPar = {'0x2005', '0x1018'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_hardcopy_protocol') ,'static');
                    vCurPar = {'0x2005', '0x1399'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_iview_bold_enabled') ,'static');
                    vCurPar = {'0x2005', '0x1019'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_inverse_reconstructed') ,'static');
                    vCurPar = {'0x2005', '0x10A2'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_is_coca') ,'static');
                    vCurPar = {'0x2005', '0x143B'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_is_coil_survey') ,'static');
                    vCurPar = {'0x2005', '0x101A'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_label_syntax') ,'static');
                    vCurPar = {'0x2005', '0x101C'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_mtc') ,'static');
                    vCurPar = {'0x2005', '0x101B'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_mp') ,'static');
                    vCurPar = {'0x2005', '0x1398'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_mobiview_enabled') ,'static');
                    vCurPar = {'0x2005', '0x1370'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_mixes_spectro') ,'static');
                    vCurPar = {'0x2005', '0x1201'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_no_film_consums') ,'static');
                    vCurPar = {'0x2005', '0x1450'}; D.set_par(vCurPar, D.get_db_value('RFR_PATIENT_no_patient_other_ids') ,'static');
                    vCurPar = {'0x2005', '0x1213'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_no_codes') ,'static');
                    vCurPar = {'0x2005', '0x1086'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_geoms') ,'static');
                    vCurPar = {'0x2005', '0x1021'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_mixes') ,'static');
                    vCurPar = {'0x2005', '0x1382'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_no_procedure_codes') ,'static');
                    vCurPar = {'0x2005', '0x1022'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_references') ,'static');
                    vCurPar = {'0x2005', '0x1199'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_no_request_excerpts') ,'static');
                    vCurPar = {'0x2005', '0x1200'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_no_sop_common') ,'static');
                    vCurPar = {'0x2005', '0x1403'}; D.set_par(vCurPar, D.get_db_value('RFR_EXAM_no_sps_codes') ,'static');
                    vCurPar = {'0x2005', '0x1023'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_slabs') ,'static');
                    vCurPar = {'0x2005', '0x1401'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_no_refd_study_sequence') ,'static');
                    vCurPar = {'0x2005', '0x1025'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_volumes') ,'static');
                    vCurPar = {'0x2005', '0x1026'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_osp') ,'static');
                    vCurPar = {'0x2005', '0x102E'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_ppg_ppu_gating') ,'static');
                    vCurPar = {'0x2005', '0x1027'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_package_mode') ,'static');
                    vCurPar = {'0x2005', '0x1029'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_partial_fourier_phase') ,'static');
                    vCurPar = {'0x2005', '0x102A'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_patient_reference_id') ,'static');
                    vCurPar = {'0x2005', '0x102B'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_percent_scan_complete') ,'static');
                    vCurPar = {'0x2005', '0x102D'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_planscan_no_images') ,'static');
                    vCurPar = {'0x2005', '0x1031'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_respiratory_gating') ,'static');
                    vCurPar = {'0x2005', '0x1381'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_scanogram_no_images') ,'static');
                    vCurPar = {'0x2005', '0x1034'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_sk') ,'static');
                    vCurPar = {'0x2001', '0x107B'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_acquisition_no') ,'static');
                    vCurPar = {'0x2005', '0x106F'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_mr_acquisition_type') ,'static');
                    vCurPar = {'0x2005', '0x142E'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_dbdt') ,'static');
                    vCurPar = {'0x2005', '0x1035'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_series_data_type') ,'static');
                    vCurPar = {'0x2005', '0x1013'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_development_mode') ,'static');
                    vCurPar = {'0x2005', '0x1014'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_diffusion') ,'static');
                    vCurPar = {'0x2001', '0x1011'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_diffusion_echo_time') ,'static');
                    vCurPar = {'0x2001', '0x1012'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_dynamic_series') ,'static');
                    vCurPar = {'0x2001', '0x1013'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_epi_factor') ,'static');
                    vCurPar = {'0x2001', '0x1025'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_echo_time_display') ,'static');
                    vCurPar = {'0x2001', '0x1082'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_echo_train_length_piim') ,'static');
                    vCurPar = {'0x2001', '0x1023'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_flip_angle') ,'static');
                    vCurPar = {'0x2005', '0x10A9'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_geometry_correction') ,'static');
                    vCurPar = {'0x2005', '0x1070'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_hardcopy_protocol_ev') ,'static');
                    vCurPar = {'0x2001', '0x1083'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_resonance_frequency') ,'static');
                    vCurPar = {'0x2001', '0x1084'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_inversion_time_piim') ,'static');
                    vCurPar = {'0x2005', '0x1036'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_series_is_cardiac') ,'static');
                    vCurPar = {'0x2001', '0x1024'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_series_is_interactive') ,'static');
                    vCurPar = {'0x2005', '0x1037'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_series_is_spectro') ,'static');
                    vCurPar = {'0x2005', '0x1431'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_local_sar') ,'static');
                    vCurPar = {'0x2001', '0x1085'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_magnetic_field') ,'static');
                    vCurPar = {'0x2005', '0x1430'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_non_proton_sar') ,'static');
                    vCurPar = {'0x2005', '0x1020'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_chemical_shifts') ,'static');
                    vCurPar = {'0x2005', '0x1414'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_diff_b_values') ,'static');
                    vCurPar = {'0x2005', '0x1415'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_diff_grad_orients') ,'static');
                    vCurPar = {'0x2001', '0x1081'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_dynamic_scans') ,'static');
                    vCurPar = {'0x2001', '0x1014'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_echoes') ,'static');
                    vCurPar = {'0x2005', '0x1428'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_label_types') ,'static');
                    vCurPar = {'0x2001', '0x1015'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_locations') ,'static');
                    vCurPar = {'0x2001', '0x1016'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_pc_directions') ,'static');
                    vCurPar = {'0x2001', '0x1017'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_phases') ,'static');
                    vCurPar = {'0x2001', '0x1018'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_slices') ,'static');
                    vCurPar = {'0x2001', '0x1060'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_stacks') ,'static');
                    vCurPar = {'0x2001', '0x1087'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_imaged_nucleus') ,'static');
                    vCurPar = {'0x2001', '0x1088'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_no_averages') ,'static');
                    vCurPar = {'0x2001', '0x1019'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_partial_matrix_scanned') ,'static');
                    vCurPar = {'0x2005', '0x1416'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_plan_mode') ,'static');
                    vCurPar = {'0x2001', '0x101B'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_prepulse_delay') ,'static');
                    vCurPar = {'0x2001', '0x101C'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_prepulse_type') ,'static');
                    vCurPar = {'0x2005', '0x142F'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_proton_sar') ,'static');
                    vCurPar = {'0x2001', '0x101D'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_reconstruction_no') ,'static');
                    vCurPar = {'0x2001', '0x101E'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_reformat_accuracy') ,'static');
                    vCurPar = {'0x2001', '0x101F'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_respiration_sync') ,'static');
                    vCurPar = {'0x2005', '0x1432'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_safety_override_mode') ,'static');
                    vCurPar = {'0x2005', '0x1033'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_scan_duration') ,'static');
                    vCurPar = {'0x2001', '0x1020'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_scanning_technique') ,'static');
                    vCurPar = {'0x2001', '0x1021'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_spir') ,'static');
                    vCurPar = {'0x2005', '0x109F'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_spectr_sel_excitation') ,'static');
                    vCurPar = {'0x2005', '0x103C'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_tone') ,'static');
                    vCurPar = {'0x2001', '0x108B'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_transmitting_coil') ,'static');
                    vCurPar = {'0x2001', '0x1022'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_water_fat_shift') ,'static');
                    vCurPar = {'0x2005', '0x102F'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_presaturation') ,'static');
                    vCurPar = {'0x2005', '0x1356'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_spectro_si_mode') ,'static');
                    vCurPar = {'0x2005', '0x1038'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_sp') ,'static');
                    vCurPar = {'0x2005', '0x1393'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_station_no') ,'static');
                    vCurPar = {'0x2005', '0x1039'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_ss') ,'static');
                    vCurPar = {'0x0032', '0x4000'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_study_comments') ,'static');
                    vCurPar = {'0x2005', '0x105F'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_study_origin') ,'static');
                    vCurPar = {'0x2005', '0x1060'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_study_sequence_no') ,'static');
                    vCurPar = {'0x2005', '0x103A'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_sub_anatomy') ,'static');
                    vCurPar = {'0x2005', '0x103B'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_trss') ,'static');
                    vCurPar = {'0x2005', '0x1400'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_volumeview_enabled') ,'static');
                    vCurPar = {'0x0018', '0x9020'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_magnetization_transfer') ,'static');
                    vCurPar = {'0x0018', '0x9011'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_multiple_spin_echo') ,'static');
                    vCurPar = {'0x0018', '0x9029'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_oversampling_phase') ,'static');
                    vCurPar = {'0x0018', '0x9078'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_parallel_acq_technique') ,'static');
                    vCurPar = {'0x0018', '0x9077'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_parallel_acquisition') ,'static');
                    vCurPar = {'0x0018', '0x9081'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_partial_fourier') ,'static');
                    vCurPar = {'0x0010', '0x0030'}; D.set_par(vCurPar, D.get_db_value('RFR_PATIENT_patient_birth_date', true) ,'static');
                    vCurPar = {'0x0010', '0x4000'}; D.set_par(vCurPar, D.get_db_value('RFR_PATIENT_patient_comments') ,'static');
                    vCurPar = {'0x0010', '0x2160'}; D.set_par(vCurPar, D.get_db_value('RFR_PATIENT_patient_ethnic_group') ,'static');
                    vCurPar = {'0x0010', '0x0020'}; D.set_par(vCurPar, D.get_db_value('RFR_PATIENT_patient_id') ,'static');
                    vCurPar = {'0x0010', '0x0010'}; D.set_par(vCurPar, D.get_db_value('RFR_PATIENT_patient_name') ,'static');
                    vCurPar = {'0x0010', '0x0040'}; D.set_par(vCurPar, D.get_db_value('RFR_PATIENT_patient_sex') ,'static');
                    vCurPar = {'0x0018', '0x9014'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_phase_contrast') ,'static');
                    vCurPar = {'0x0018', '0x0095'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_pixel_bandwidth') ,'static');
                    vCurPar = {'0x0008', '0x9205'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_pixel_presentation') ,'static');
                    vCurPar = {'0x2005', '0x1447'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_power_optimization') ,'static');
                    vCurPar = {'0x1001', '0x1007'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_print_status') ,'static');
                    vCurPar = {'0x0018', '0x9005'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_pulse_sequence_name') ,'static');
                    vCurPar = {'0x0018', '0x9240'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_rf_echo_train_length') ,'static');
                    vCurPar = {'0x1001', '0x100C'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_remote_activ_status') ,'static');
                    vCurPar = {'0x0018', '0x9170'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_resp_motion_comp_techn') ,'static');
                    vCurPar = {'0x0020', '0x9254'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_resp_interval_time') ,'static');
                    vCurPar = {'0x2005', '0x1032'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_sample_representation') ,'static');
                    vCurPar = {'0x0018', '0x9024'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_saturation_recovery') ,'static');
                    vCurPar = {'0x2001', '0x1062'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_committed') ,'static');
                    vCurPar = {'0x2001', '0x10CC'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_derivation_desc') ,'static');
                    vCurPar = {'0x0018', '0x5100'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_patient_position') ,'static');
                    vCurPar = {'0x2005', '0x10C0'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_series_scan_sequence') ,'static');
                    vCurPar = {'0x2005', '0x1134'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_transaction_uid') ,'static');
                    vCurPar = {'0x2001', '0x1061'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_transmitted') ,'static');
                    vCurPar = {'0x2001', '0x106E'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_series_type') ,'static');
                    vCurPar = {'0x2005', '0x101E'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_mip_protocol') ,'static');
                    vCurPar = {'0x2005', '0x101F'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_mpr_protocol') ,'static');
                    vCurPar = {'0x0018', '0x9027'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_spatial_pre_saturation') ,'static');
                    vCurPar = {'0x0018', '0x9016'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_spoiling') ,'static');
                    vCurPar = {'0x0018', '0x9017'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_steady_state_pulse_seq') ,'static');
                    vCurPar = {'0x1001', '0x100A'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_storage_commit_status') ,'static');
                    vCurPar = {'0x0008', '0x0020'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_study_date', true) ,'static');
                    vCurPar = {'0x0008', '0x1030'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_study_description') ,'static');
                    vCurPar = {'0x0020', '0x0010'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_study_id_dicom', true) ,'static');
                    vCurPar = {'0x0020', '0x000D'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_study_instance_uid') ,'static');
                    vCurPar = {'0x0038', '0x0050'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_patient_special_needs') ,'static');
                    vCurPar = {'0x0038', '0x0500'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_patient_state') ,'static');
                    vCurPar = {'0x0010', '0x2180'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_patient_occupation') ,'static');
                    vCurPar = {'0x0010', '0x1020'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_study_patients_size') ,'static');
                    vCurPar = {'0x0010', '0x1030'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_patient_weight') ,'static');
                    vCurPar = {'0x0008', '0x0030'}; D.set_par(vCurPar, D.get_db_value('RFR_STUDY_study_time', true) ,'static');
                    vCurPar = {'0x0018', '0x9021'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_t2_preparation') ,'static');
                    vCurPar = {'0x2005', '0x1444'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_tfe_factor') ,'static');
                    vCurPar = {'0x0018', '0x9019'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_tag_angle_first_axis') ,'static');
                    vCurPar = {'0x0018', '0x9030'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_tag_spacing_first_dim') ,'static');
                    vCurPar = {'0x0018', '0x9028'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_tagging') ,'static');
                    vCurPar = {'0x0018', '0x9051'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_transmit_coil_type') ,'static');
                    vCurPar = {'0x2005', '0x1426'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_viewing_hardcopy_only') ,'static');
                    vCurPar = {'0x0008', '0x9206'}; D.set_par(vCurPar, D.get_db_value('RFR_SERIES_volumetric_properties') ,'static');
                end
                
            end
        end
        function dicom_get_dynamic_values_from_parfile(D, Parfile, image_nr)
            
            sop_uid = D.get_uid('');
            vCurPar = 'SOPInstanceUID';  D.set_par(vCurPar, sop_uid,'dynamic');
            vCurPar = 'MediaStorageSOPInstanceUID';  D.set_par(vCurPar, sop_uid,'dynamic');
            
            vCurPar = 'SliceLocation';  D.set_par(vCurPar, Parfile.ImageInformation.SliceNumber(image_nr));
            vCurPar = {'0x2001', '0x100a'};  D.set_par(vCurPar,Parfile.ImageInformation.SliceNumber(image_nr));
            vCurPar = 'EchoNumbers';  D.set_par(vCurPar, Parfile.ImageInformation.EchoNumber(image_nr));
            vCurPar = 'TemporalPositionIdentifier';  D.set_par(vCurPar, Parfile.ImageInformation.DynamicScanNumber(image_nr));
            vCurPar = {'0x2001', '0x1008'};  D.set_par(vCurPar, Parfile.ImageInformation.CardiacPhaseNumber(image_nr));
            
            switch (Parfile.ImageInformation.ImageTypeMr(image_nr))
                case 0
                    cur_type_str = 'M';
                case 1
                    cur_type_str = 'R';
                case 2
                    cur_type_str = 'I';
                case 3
                    cur_type_str = 'P';
            end
            
            image_type = ['ORIGINAL\PRIMARY\',  cur_type_str, '_', strtrim(Parfile.Technique), '\', cur_type_str, '\', strtrim(Parfile.Technique)];
            
            vCurPar = {'0x2005', '0x1011'};  D.set_par(vCurPar, cur_type_str,'dynamic');
            vCurPar = 'ImageType';  D.set_par(vCurPar, image_type,'dynamic');
            if( Parfile.ImageInformation.ScanningSequence(image_nr) == 0 )
                vCurPar = {'0x2005', '0x106e'};  D.set_par(vCurPar, 'FFE', 'dynamic');	  % scanning sequence
            else
                vCurPar = {'0x2005', '0x106e'};  D.set_par(vCurPar, 'SE', 'dynamic');	  % scanning sequence
            end
            
            vCurPar = 'InstanceNumber';  D.set_par(vCurPar,Parfile.ImageInformation.IndexInRECFile(image_nr) + 1);
            vCurPar = 'BitsAllocated';  D.set_par(vCurPar,Parfile.ImageInformation.ImagePixelSize(image_nr) );
            vCurPar = 'PercentSampling';  D.set_par(vCurPar,Parfile.ImageInformation.ScanPercentage(image_nr));
            vCurPar = 'Rows';  D.set_par(vCurPar,Parfile.ImageInformation.ReconResolution(image_nr,1));
            vCurPar = 'Columns';  D.set_par(vCurPar,Parfile.ImageInformation.ReconResolution(image_nr,2));
%             vCurPar = 'RescaleIntercept';  D.set_par(vCurPar,Parfile.ImageInformation.RescaleIntercept(image_nr));
%             vCurPar = 'RescaleSlope';  D.set_par(vCurPar,Parfile.ImageInformation.RescaleSlope(image_nr));
%             vCurPar = 'RescaleType';  D.set_par(vCurPar, 'normalized','dynamic');
            vCurPar = 'RescaleIntercept';  D.set_par(vCurPar,Parfile.ImageInformation.RescaleIntercept(image_nr)*1e-3*max(Parfile.PhaseEncodingVelocity)/pi); % LMG edit: add dicom tag for circle
            vCurPar = 'RescaleSlope';  D.set_par(vCurPar,Parfile.ImageInformation.RescaleSlope(image_nr)*1e-3*max(Parfile.PhaseEncodingVelocity)/pi); % LMG edit: add dicom tag for circle
            vCurPar = 'RescaleType';  D.set_par(vCurPar, 'cm/s','dynamic'); % LMG edit: add dicom tag for circle
            vCurPar = 'RescaleInterceptOriginal';  D.set_par(vCurPar,Parfile.ImageInformation.RescaleIntercept(image_nr)*1e-3*max(Parfile.PhaseEncodingVelocity)/pi); % LMG edit: add dicom tag for circle
            vCurPar = 'RescaleSlopeOriginal';  D.set_par(vCurPar,Parfile.ImageInformation.RescaleSlope(image_nr)*1e-3*max(Parfile.PhaseEncodingVelocity)/pi); % LMG edit: add dicom tag for circle
            vCurPar = 'RescaleTypeOriginal';  D.set_par(vCurPar, 'cm/s','dynamic'); % LMG edit: add dicom tag for circle
            vCurPar = 'RealWorldValueMappingSequence';  D.set_par(vCurPar, ''); % LMG edit: add dicom tag for circle
            vCurPar = 'LUTLabel';  D.set_par(vCurPar, 'Philips'); % LMG edit: add dicom tag for circle
            vCurPar = 'RealWorldValueLastValueMapped';  D.set_par(vCurPar,4096); % LMG edit: add dicom tag for circle
            vCurPar = 'RealWorldValueFirstValueMapped';  D.set_par(vCurPar,0); % LMG edit: add dicom tag for circle
            vCurPar = 'RealWorldValueSlope';  D.set_par(vCurPar,Parfile.ImageInformation.RescaleSlope(image_nr)*1e-3*max(Parfile.PhaseEncodingVelocity)/pi); % LMG edit: add dicom tag for circle
            vCurPar = 'RealWorldValueIntercept';  D.set_par(vCurPar,Parfile.ImageInformation.RescaleIntercept(image_nr)*1e-3*max(Parfile.PhaseEncodingVelocity)/pi); % LMG edit: add dicom tag for circle
            vCurPar = 'WindowCenter';  D.set_par(vCurPar,Parfile.ImageInformation.WindowCenter(image_nr));
            vCurPar = 'WindowWidth';  D.set_par(vCurPar,Parfile.ImageInformation.WindowWidth(image_nr));
            if Parfile.ImageInformation.ImageTypeMr(image_nr) == 0
                slope = Parfile.ImageInformation.RescaleSlope(image_nr)*1e-3*max(Parfile.PhaseEncodingVelocity)/pi;
                intercept = Parfile.ImageInformation.RescaleIntercept(image_nr)*1e-3*max(Parfile.PhaseEncodingVelocity)/pi;
                vCurPar = 'WindowCenter';  D.set_par(vCurPar,2048*slope-intercept); % LMG edit: add dicom tag for circle
                vCurPar = 'WindowWidth';  D.set_par(vCurPar,4096*slope-intercept); % LMG edit: add dicom tag for circle
            elseif Parfile.ImageInformation.ImageTypeMr(image_nr) == 3
                vCurPar = 'WindowCenter';  D.set_par(vCurPar,0); % LMG edit: add dicom tag for circle
                vCurPar = 'WindowWidth';  D.set_par(vCurPar,max(Parfile.PhaseEncodingVelocity)*2); % LMG edit: add dicom tag for circle
            end
            vCurPar = {'0x2005', '0x1000'};  D.set_par(vCurPar,Parfile.ImageInformation.ImageAngulation(image_nr, 1));
            vCurPar = {'0x2005', '0x1001'};  D.set_par(vCurPar,Parfile.ImageInformation.ImageAngulation(image_nr, 2));
            vCurPar = {'0x2005', '0x1002'};  D.set_par(vCurPar,Parfile.ImageInformation.ImageAngulation(image_nr, 3));
            vCurPar = {'0x2005', '0x1008'};  D.set_par(vCurPar,Parfile.ImageInformation.ImageOffcentre(image_nr, 1));
            vCurPar = {'0x2005', '0x1009'};  D.set_par(vCurPar,Parfile.ImageInformation.ImageOffcentre(image_nr, 2));
            vCurPar = {'0x2005', '0x100a'};  D.set_par(vCurPar,Parfile.ImageInformation.ImageOffcentre(image_nr, 3));
            
            [aDirectionCosinesRL, aDirectionCosinesAP, aDirectionCosinesFH] = DICOMExporter.get_direction_cosines(Parfile, image_nr);
            direction_cosine = [aDirectionCosinesRL(1); aDirectionCosinesAP(1); aDirectionCosinesFH(1); aDirectionCosinesRL(2); aDirectionCosinesAP(2); aDirectionCosinesFH(2)];
            
            vCurPar = 'ImageOrientationPatient';  D.set_par(vCurPar, direction_cosine,'dynamic');
            
            aDirectionCosinesX = [aDirectionCosinesRL(1); aDirectionCosinesAP(1); aDirectionCosinesFH(1)];
            aDirectionCosinesY = [aDirectionCosinesRL(2); aDirectionCosinesAP(2); aDirectionCosinesFH(2)];
            aDirectionCosinesZ = [aDirectionCosinesRL(3); aDirectionCosinesAP(3); aDirectionCosinesFH(3)];
            
            center_to_upper_left_corner_RAF = zeros(3, 1);
            center_to_upper_left_corner_RAF(1) = -(Parfile.ImageInformation.PixelSpacing(image_nr,1)*Parfile.ImageInformation.ReconResolution(image_nr, 1) / 2.0) * aDirectionCosinesX(1);
            center_to_upper_left_corner_RAF(2) = -(Parfile.ImageInformation.PixelSpacing(image_nr,1)*Parfile.ImageInformation.ReconResolution(image_nr, 1) / 2.0) * aDirectionCosinesX(2);
            center_to_upper_left_corner_RAF(3) = -(Parfile.ImageInformation.PixelSpacing(image_nr,1)*Parfile.ImageInformation.ReconResolution(image_nr, 1) / 2.0) * aDirectionCosinesX(3);
            
            center_to_upper_left_corner_RAF(1) = center_to_upper_left_corner_RAF(1) -(Parfile.ImageInformation.PixelSpacing(image_nr, 2)*Parfile.ImageInformation.ReconResolution(image_nr, 2) / 2.0) * aDirectionCosinesY(1);
            center_to_upper_left_corner_RAF(2) = center_to_upper_left_corner_RAF(2) -(Parfile.ImageInformation.PixelSpacing(image_nr, 2)*Parfile.ImageInformation.ReconResolution(image_nr, 2) / 2.0) * aDirectionCosinesY(2);
            center_to_upper_left_corner_RAF(3) = center_to_upper_left_corner_RAF(3) -(Parfile.ImageInformation.PixelSpacing(image_nr, 2)*Parfile.ImageInformation.ReconResolution(image_nr, 2) / 2.0) * aDirectionCosinesY(3);
            
            upper_left_corner(1) = Parfile.ImageInformation.ImageOffcentre(image_nr, 3) + center_to_upper_left_corner_RAF(1);
            upper_left_corner(2) = Parfile.ImageInformation.ImageOffcentre(image_nr, 1) + center_to_upper_left_corner_RAF(2);
            upper_left_corner(3) = Parfile.ImageInformation.ImageOffcentre(image_nr, 2) + center_to_upper_left_corner_RAF(3);
            
            image_position_patient = [upper_left_corner(1); upper_left_corner(2); upper_left_corner(3)];
            vCurPar = 'ImagePositionPatient';  D.set_par(vCurPar, image_position_patient,'dynamic');
            vCurPar = 'SliceThickness';  D.set_par(vCurPar,Parfile.ImageInformation.SliceThickness(image_nr));
            vCurPar = 'SpacingBetweenSlices';  D.set_par(vCurPar, Parfile.ImageInformation.SliceThickness(image_nr) + Parfile.ImageInformation.SliceGap(image_nr));
            if( Parfile.ImageInformation.SliceOrientation(image_nr) == 1 )
                vCurPar = {'0x2001', '0x100b'};  D.set_par(vCurPar, 'TRANSVERSAL' );
            elseif( Parfile.ImageInformation.SliceOrientation(image_nr) == 2 )
                vCurPar = {'0x2001', '0x100b'};  D.set_par(vCurPar, 'SAGITTAL' );
            elseif( Parfile.ImageInformation.SliceOrientation(image_nr) == 3 )
                vCurPar = {'0x2001', '0x100b'};  D.set_par(vCurPar, 'CORONAL' );
            end
            
            vCurPar = {'0x2005', '0x1063'};  D.set_par(vCurPar,Parfile.ImageInformation.FmriStatusIndication(image_nr)); 	  % fMRI status indication
            
            vPixelSpacing = [Parfile.ImageInformation.PixelSpacing(image_nr, 1); Parfile.ImageInformation.PixelSpacing(image_nr, 2)];
            vCurPar = 'PixelSpacing';  D.set_par(vCurPar, vPixelSpacing,'dynamic');
            vCurPar = 'ImagerPixelSpacing';  D.set_par(vCurPar, vPixelSpacing,'dynamic');
            vCurPar = 'EchoTime';  D.set_par(vCurPar,Parfile.ImageInformation.EchoTime(image_nr));
            
            acq_time_str = '120000.65000';
            date_time = D.format_date_time(Parfile.ExaminationDateTime);
            if (~isempty(date_time(2)))
                acq_time_str = date_time(2);
            end
            
            dyn_begin_time = Parfile.ImageInformation.DynScanBeginTime(image_nr) - Parfile.ImageInformation.DynScanBeginTime(1);
            
            try
                acq_time_final_str = datestr(datetime(acq_time_str, 'InputFormat', 'HHmmss.SSSSS') + seconds(dyn_begin_time), 'HHMMSS.FFF00');
            catch
                acq_time_final_str = acq_time_str;
            end
            
            dyn_scan_time = 0;
            if (Parfile.DynamicScan == 1)
                dyn_scan_time = Parfile.ImageInformation.DynScanBeginTime(image_nr);
                vCurPar = 'TriggerTime';  D.set_par(vCurPar,1000*dyn_scan_time);
            else
                vCurPar = 'TriggerTime';  D.set_par(vCurPar, Parfile.ImageInformation.TriggerTime(image_nr));
            end
            vCurPar = {'0x2005', '0x10A0'};  D.set_par(vCurPar, dyn_scan_time / 1000.0 );
            vCurPar = 'AcquisitionTime';  D.set_par(vCurPar, acq_time_final_str,'dynamic');
            vCurPar = 'ContentTime';  D.set_par(vCurPar, acq_time_final_str,'dynamic');
            
            % get the diffusion parameter
            diffusion_b_factor = 0;
            diffusion_b_value_number = 1;
            gradient_orientation_number = 1;
            diffusionAP = 0;
            diffusionFH = 0;
            diffusionRL = 0;
            diffusionM = 0;
            diffusionP = 0;
            diffusionS = 0;
            DiffusionBFactor = 0;
            DiffusionDirection = '';
            
            if (Parfile.Diffusion == 1)
                diffusion_b_value_number = Parfile.ImageInformation.DiffusionBValueNumber(image_nr);
                gradient_orientation_number = Parfile.ImageInformation.GradientOrientationNumber(image_nr);
                
                diffusion_b_factor = Parfile.ImageInformation.DiffusionBFactor(image_nr);
                DiffusionBFactor = diffusion_b_factor;
                diffusionAP = Parfile.ImageInformation.Diffusion(image_nr, 1);
                diffusionFH = Parfile.ImageInformation.Diffusion(image_nr, 2);
                diffusionRL = Parfile.ImageInformation.Diffusion(image_nr, 3);
                
                switch (Parfile.ImageInformation.get_SliceOrientation(image_nr))
                    case 3 %COR:
                        diffusionM = diffusionFH;  diffusionP = diffusionRL; diffusionS = diffusionAP;
                    case 2 %SAG:
                        diffusionM = diffusionFH;  diffusionP = diffusionAP; diffusionS = diffusionRL;
                    case 1 %TRA:
                        diffusionM = diffusionAP;  diffusionP = diffusionRL; diffusionS = diffusionFH;
                end
                
                DiffusionMPS = [diffusionM; diffusionP; diffusionS];
                
                DiffusionDirection = 'O';
                if (diffusionM == 0 && diffusionP == 0 && abs(diffusionS) == 1)
                    DiffusionDirection = 'S';
                elseif (diffusionM == 0 && abs(diffusionP) == 1 && diffusionS == 0)
                    DiffusionDirection = 'P';
                elseif (abs(diffusionM) == 1 && diffusionP == 0 && diffusionS == 0)
                    DiffusionDirection = 'M';
                elseif (diffusionM == 0 && diffusionP == 0 && diffusionS == 0)
                    DiffusionDirection = 'I';
                end
            else
                DiffusionMPS = [0; 0; 0];
                DiffusionBFactor = 0;
            end
            
            vCurPar = {'0x2005', '0x10b0'};  D.set_par(vCurPar,diffusionRL);
            vCurPar = {'0x2005', '0x10b1'};  D.set_par(vCurPar,diffusionAP);
            vCurPar = {'0x2005', '0x10b2'};  D.set_par(vCurPar,diffusionFH);
            vCurPar = 'DiffusionBValue';  D.set_par(vCurPar,diffusion_b_factor);
            vCurPar = 'DiffusionGradientOrientation';  D.set_par(vCurPar, DiffusionMPS,'dynamic');
            vCurPar = {'0x2001', '0x1003'};  D.set_par(vCurPar,DiffusionBFactor);
            vCurPar = {'0x2001', '0x1004'};  D.set_par(vCurPar, DiffusionDirection,'dynamic');
            vCurPar = 'NumberOfAverages';  D.set_par(vCurPar,Parfile.ImageInformation.NumberOfAverages(image_nr));
            vCurPar = 'FlipAngle';  D.set_par(vCurPar,Parfile.ImageInformation.ImageFlipAngle(image_nr));
            
            if Parfile.ImageInformation.CardiacFrequency(image_nr) > 0
                vCurPar = 'HeartRate';  D.set_par(vCurPar,Parfile.ImageInformation.CardiacFrequency(image_nr));
            end
            
            if Parfile.ImageInformation.MinRRInterval(image_nr) > 0
                vCurPar = 'LowRRValue';  D.set_par(vCurPar,Parfile.ImageInformation.MinRRInterval(image_nr));
            end
            if Parfile.ImageInformation.MaxRRInterval(image_nr) > 0
                vCurPar = 'HighRRValue';  D.set_par(vCurPar,Parfile.ImageInformation.MaxRRInterval(image_nr));
            end
            vCurPar = {'0x2005', '0x1444'};  D.set_par(vCurPar,Parfile.ImageInformation.TurboFactor(image_nr));
            vCurPar = 'InversionTime';  D.set_par(vCurPar,Parfile.ImageInformation.InversionDelay(image_nr));
            vCurPar = 'PhotometricInterpretation';  D.set_par(vCurPar, 'MONOCHROME2','dynamic');
            
            vBitsStored = 12;
            vPixelRepresentation = 0;
            vCurPar = 'BitsStored';  D.set_par(vCurPar, vBitsStored);
            vCurPar = 'HighBit';  D.set_par(vCurPar, vBitsStored - 1);
            vCurPar = 'PixelRepresentation';  D.set_par(vCurPar,vPixelRepresentation);
            vCurPar = {'0x2001', '0x1006'};  D.set_par(vCurPar, 'N','dynamic');
        end
        function dicom_get_dynamic_values(D, Parfile, image_nr, cur_type, file_number)
            D.dicom_get_dynamic_values_from_parfile(Parfile, image_nr );
            
            try
                % get series number (changes with the flow direction)
                series_number = str2double(D.mGoalcPars.GetValue('RFR_SERIES_DICOM_SERIES_NUMBER')) + file_number;
                vCurPar = {'0x0020', '0x0011'};  D.set_par(vCurPar,series_number);
                
                vCurPar = 'MRAcquisitionPhaseEncodingStepsOutOfPlane';  D.set_par(vCurPar, str2double(D.mGoalcPars.GetValue('RFR_SERIES_PIIM_MR_SERIES_NR_OF_SLICES')));
                
                if (strcmpi(D.mGoalcPars.GetValue('RFR_SERIES_PIIM_SERIES_SCAN_SEQUENCE'),'GR'))
                    scanning_sequence = {'FFE', 'SE'};
                else
                    scanning_sequence = {'SE', 'FFE'};
                end
                
                if (~isempty(strfind(D.mGoalcPars.GetValue('RFR_SERIES_PIIM_MR_SERIES_SCANNING_TECHNIQUE_DESC'), 'MIX')) && ...
                        Parfile.ImageInformation.ScanningSequence(image_nr) < length(scanning_sequence) )
                    vCurPar = {'0x2005', '0x106e'};  D.set_par(vCurPar,scanning_sequence{Parfile.ImageInformation.ScanningSequence(image_nr)}); 	  % scanning sequence
                else
                    vCurPar = {'0x2005', '0x106e'};  D.set_par(vCurPar,scanning_sequence{1}); 	  % scanning sequence
                end
                
                % if we have the database parameter then use the time from there
                acq_time_str = D.mGoalcPars.GetValue('RFR_SERIES_DICOM_SERIES_TIME');
                dyn_scan_time = 0.0;
                if (Parfile.DynamicScan == 1)
                    dyn_scan_time = Parfile.ImageInformation.DynScanBeginTime(image_nr);
                end
                acq_time_final_str = datestr(datetime(acq_time_str, 'InputFormat', 'HHmmss.SSSSS') + seconds(dyn_scan_time), 'HHMMSS.FFF00');
                
                vCurPar = 'AcquisitionTime';  D.set_par(vCurPar,acq_time_final_str);
                vCurPar = 'ContentTime';  D.set_par(vCurPar,acq_time_final_str);
                
                if strcmp(D.mRelease, 'REL53')
                    name = 'VAL01_ACQ_echoes';
                else
                    name = 'UGN1_ACQ_echoes';
                end
                if (D.mGoalcPars.IsParameter(name))
                    
                    % Oversampling in Phase direction (seemes strange the way it is but its the same as Philips does it)
                    if strcmp(D.mRelease, 'REL53')
                        name = 'VAL01_ACQ_echoes';
                    else
                        name = 'UGN1_ACQ_echoes';
                    end
                    nr_echoes = D.mGoalcPars.GetValue(name);
                    if strcmp(D.mRelease, 'REL53')
                        name = 'VAL01_ACQ_mixes';
                    else
                        name = 'UGN1_ACQ_mixes';
                    end
                    nr_mixes = D.mGoalcPars.GetValue(name);
                    kx_ovs = D.mGoalcPars.GetValue('RC_oversample_factors', 1);
                    ky_ovs = D.mGoalcPars.GetValue('RC_oversample_factors', nr_echoes*nr_mixes);
                    kz_ovs = D.mGoalcPars.GetValue('RC_oversample_factors', 2 * nr_echoes*nr_mixes);
                    OversamplingPhase = 'NONE';
                    if (ky_ovs > 1)
                        OversamplingPhase = '2D';
                        if (kz_ovs > 1)
                            OversamplingPhase = '3D';
                        end
                    end
                    vCurPar = 'OversamplingPhase';  D.set_par(vCurPar,OversamplingPhase);
                    
                    
                    cur_type_str = '';
                    switch (Parfile.ImageInformation.ImageTypeMr(image_nr))
                        case 0
                            cur_type_str = 'M';
                        case 1
                            cur_type_str = 'R';
                        case 2
                            cur_type_str = 'I';
                        case 3
                            cur_type_str = 'P';
                    end
                    
                    image_type = ['ORIGINAL\PRIMARY\', cur_type_str, '_', D.mGoalcPars.GetValue('EX_ACQ_imaging_sequence'),  '\', cur_type_str, '\', D.mGoalcPars.GetValue('EX_ACQ_imaging_sequence')];
                    % LMG start: correct phase contrast image type
                    if ( D.mGoalcPars.GetValue('EX_PC_angio_mode', [], 1) == 2 ) && strcmp(image_type,'ORIGINAL\PRIMARY\P_FFE\P\FFE')
                        image_type = 'ORIGINAL\PRIMARY\PHASE CONTRAST M\P\PCA';
                    end
                    % LMG end
                    vCurPar = 'ImageType';  D.set_par(vCurPar,image_type);
                    
                    
                    % find out the current location (somehow) and then get the in plane phase encoding direction (as a workaround we only get the first value)
                    if strcmp(D.mRelease, 'REL53')
                        name = 'VAL01_ACQ_scan_mode';
                    else
                        name = 'UGN1_ACQ_scan_mode';
                    end
                    is3d = D.mGoalcPars.GetValue(name, [], 1) == 1;
                    cur_loca = Parfile.ImageInformation.SliceNumber(image_nr);
                    if (is3d)
                        cur_loca = 1;
                    end
                    
                    if( D.mGoalcPars.GetValue('RC_phase_encoding_directions', cur_loca, 1))
                        vCurPar = 'InPlanePhaseEncodingDirection';  D.set_par(vCurPar, 'COL');
                    else
                        vCurPar = 'InPlanePhaseEncodingDirection';  D.set_par(vCurPar, 'ROW');
                    end
                    
                    % number of slices per stack
                    nr_stacks = D.mGoalcPars.GetValue('PS_nr_of_stacks');
                    cur_stack = floor((Parfile.ImageInformation.SliceNumber(image_nr)-1)/nr_stacks);
                    stack_slices = D.mGoalcPars.GetValue('EX_GEO_stacks_slices');
                    slices_per_stack = Parfile.MaxNumberOfSlicesLocations;
                    if (cur_stack >= 0 && cur_stack < length(stack_slices))
                        slices_per_stack = stack_slices(cur_stack+1);
                    end
                    if D.mGoalcPars.IsParameter('PS_stack_orientations')
                        try
                            stack_ori = D.mGoalcPars.GetValue('PS_stack_orientations', cur_stack+1);
                        catch
                            stack_ori = D.mGoalcPars.GetValue('PS_stack_orientations', 1);
                        end
                    end
                    vCurPar = {'0x2005', '0x1081'};  D.set_par(vCurPar,stack_ori);
                    vCurPar = {'0x2001', '0x102d'};  D.set_par(vCurPar,slices_per_stack);
                    
                    % Not correct, but I don't understand the Philips values yet
                    vCurPar = 'SliceLocation';  D.set_par(vCurPar,Parfile.ImageInformation.SliceNumber(image_nr));
                end
            catch
                warning('Could not get all DICOM tags. Please make sure you have the most recent version of the ReconFrame patch');
            end
        end
        function set_par( D, par_name, value, type )
            if( nargin == 3 )
                type = 'dynamic';
            end
            if( isnan(value ) )
                return;
            end
            if( iscell(par_name) )
                group = par_name{1};
                element = par_name{2};
                group = strrep(group, '0x', '');
                element = strrep(element, '0x', '');
                par_name = dicomlookup(group, element);
                if( isempty(par_name))
                    return;
                end
                group_dec = hex2dec(group);
                element_dec = hex2dec(element);
            else
                if( isempty(dicomlookup(par_name)))
                    return;
                end
                [group_dec, element_dec] = dicomlookup(par_name);
            end
            
            if( ischar(value) )
                value = strtrim(value);
            end
            
            % check if the current tag is in one or multiple sequences
            sequence_ind = find( D.mSequences{1}(:,1) == group_dec & D.mSequences{1}(:,2) == element_dec , 1);
            if( ~isempty(sequence_ind))
                sequences = D.mSequences{2}(sequence_ind);
                sequences = sequences{1};
                for i = 1:size(sequences, 1)
                    
                    if( all(sequences(i,:)) == 0 )
                        if( strcmpi(type, 'static' ))
                            D.mInfoStatic.(par_name) = value;
                        else
                            D.mInfo.(par_name) = value;
                        end
                    else
                        % get the sequence names
                        sequence_name = dicomlookup(dec2hex(sequences(i,1)), dec2hex(sequences(i,2)));
                        
                        if( strcmpi(type, 'static' ))
                            D.mInfoStatic.(sequence_name).Item_1.(par_name) = value;
                        else
                            D.mInfo.(sequence_name).Item_1.(par_name) = value;
                        end
                    end
                end
            else
                if( strcmpi(type, 'static' ))
                    D.mInfoStatic.(par_name) = value;
                else
                    D.mInfo.(par_name) = value;
                end
            end
        end
        function value = get_db_value( D, name, as_string)
            if( nargin == 2 )
                as_string = false;
            end
            try
                if( iscell(name) )
                    
                    if( as_string )
                        for i = 1:length(name)
                            p = D.mGoalcPars.GetParameter(name{i});
                            value_temp = p.GetValue();
                            if( i == 1)
                                value = value_temp;
                            else
                                value = [value, '\',  value_temp];
                            end
                        end
                    else
                        for i = 1:length(name)
                            p = D.mGoalcPars.GetParameter(name{i});                                                     
                            if isnumeric(p.GetValue())
                                value(i) = p.GetValue();
                            else
                                value_temp = str2double(p.GetValue());
                                if( isnan(value_temp))
                                    if( i == 1)
                                        value = p.GetValue();
                                    else
                                        value = [value, '\',  p.GetValue()];
                                    end
                                else
                                    value(i) = value_temp;
                                end
                            end
                        end
                    end
                    
                else
                    p = D.mGoalcPars.GetParameter(name);
                    if( as_string )
                        value = p.GetValue();
                    else
                        if isnumeric(p.GetValue())
                            value = p.GetValue();
                        else
                            value = str2double(p.GetValue());
                            if( isnan(value))
                                value = p.GetValue();
                            end
                        end
                    end
                end
            catch
                value = NaN;
            end
        end
        function get_filenames(D, nr_recfiles)
            if( isempty(D.mGoalcPars) )
                if (nr_recfiles > 1)
                    for  i = 1:nr_recfiles
                        D.mOutDirRelativePath{i} = num2str(i);
                    end
                else
                    D.mOutDirRelativePath{1} = '';
                end
            else
                try
                    if (strcmpi(D.mRelease, 'REL3'))
                        venc(1) = D.mGoalcPars.GetValue('RFR_SERIES_pc_velocity0');
                        venc(2) = D.mGoalcPars.GetValue('RFR_SERIES_pc_velocity1');
                        venc(3) = D.mGoalcPars.GetValue('RFR_SERIES_pc_velocity2');
                    else
                        venc(1) = str2double(D.mGoalcPars.GetValue('RFR_SERIES_PIIM_MR_SERIES_PC_VELOCITY0'));
                        venc(2) = str2double(D.mGoalcPars.GetValue('RFR_SERIES_PIIM_MR_SERIES_PC_VELOCITY1'));
                        venc(3) = str2double(D.mGoalcPars.GetValue('RFR_SERIES_PIIM_MR_SERIES_PC_VELOCITY2'));
                    end
                    
                    isFlow = any(venc > 0);
                    is3dflow = isFlow && all( venc > 0 );
                    if strcmp(D.mRelease, 'REL53')
                        ismultivenc = D.mGoalcPars.IsParameter('VAL01_PC_multi_venc') && D.mGoalcPars.GetValue('VAL01_PC_multi_venc', [], 1) == 1;
                    else
                        ismultivenc = D.mGoalcPars.IsParameter('UGN1_PC_multi_venc') && D.mGoalcPars.GetValue('UGN1_PC_multi_venc', [], 1) == 1;
                    end
                    
                    
                    if (nr_recfiles > 1)
                        if (isFlow && is3dflow && ~ismultivenc)
                            if(nr_recfiles == 1)
                                D.mOutDirRelativePath{1} = 'M';
                            elseif (nr_recfiles == 2)
                                D.mOutDirRelativePath{1} = 'M';
                                D.mOutDirRelativePath{2} = 'P';
                            elseif (nr_recfiles == 3)
                                D.mOutDirRelativePath{1} = 'M';
                                D.mOutDirRelativePath{2} = 'P';
                                D.mOutDirRelativePath{3} = 'S';
                            else
                                for i = 1:nr_recfiles
                                    D.mOutDirRelativePath{i} = num2str(i);
                                end
                            end
                        else
                            for i = 1:nr_recfiles
                                D.mOutDirRelativePath{i} = num2str(i);
                            end
                        end
                    else
                        D.mOutDirRelativePath{1} = '';
                    end
                catch
                    if (nr_recfiles > 1)
                        for  i = 1:nr_recfiles
                            D.mOutDirRelativePath{i} = num2str(i);
                        end
                    else
                        D.mOutDirRelativePath{1} = '';
                    end
                end
            end
        end
        function define_sequences(D)
            tags = [];
            sequences = [];
            loop = 1;
            
            % Define the Performed Procedure Step Sequence
            tags(loop,:) = [hex2dec('0008'), hex2dec('0014')]; sequences{loop} = [hex2dec('8'), hex2dec('1111')]; loop = loop+1;
            tags(loop,:) = [hex2dec('0008'), hex2dec('1150')]; sequences{loop} = [hex2dec('8'), hex2dec('1111')]; loop = loop+1;
            tags(loop,:) = [hex2dec('0008'), hex2dec('1155')]; sequences{loop} = [hex2dec('8'), hex2dec('1111')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('0014')]; sequences{loop} = [hex2dec('8'), hex2dec('1111')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1404')]; sequences{loop} = [hex2dec('8'), hex2dec('1111')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1406')]; sequences{loop} = [hex2dec('8'), hex2dec('1111')]; loop = loop+1;
            % Add the following tags to the base sequence as well
            tags(loop,:) = [hex2dec('0008'), hex2dec('0012')]; sequences{loop} = [hex2dec('8'), hex2dec('1111'); hex2dec('0'), hex2dec('0')]; loop = loop+1;
            tags(loop,:) = [hex2dec('0008'), hex2dec('0013')]; sequences{loop} = [hex2dec('8'), hex2dec('1111'); hex2dec('0'), hex2dec('0')]; loop = loop+1;
            tags(loop,:) = [hex2dec('0020'), hex2dec('0013')]; sequences{loop} = [hex2dec('8'), hex2dec('1111'); hex2dec('0'), hex2dec('0')]; loop = loop+1;
            
            
            % Define the ReferencedImageSequence
            tags(loop,:) = [hex2dec('0008'), hex2dec('1150')]; sequences{loop} = [hex2dec('8'), hex2dec('1140')]; loop = loop+1;
            tags(loop,:) = [hex2dec('0008'), hex2dec('1155')]; sequences{loop} = [hex2dec('8'), hex2dec('1140')]; loop = loop+1;
            
            
            % Define the PerformedProtocolCodeSequence
            tags(loop,:) = [hex2dec('0008'), hex2dec('0100')]; sequences{loop} = [hex2dec('40'), hex2dec('260')]; loop = loop+1;
            tags(loop,:) = [hex2dec('0008'), hex2dec('0102')]; sequences{loop} = [hex2dec('40'), hex2dec('260')]; loop = loop+1;
            tags(loop,:) = [hex2dec('0008'), hex2dec('0104')]; sequences{loop} = [hex2dec('40'), hex2dec('260')]; loop = loop+1;
            tags(loop,:) = [hex2dec('0008'), hex2dec('010b')]; sequences{loop} = [hex2dec('40'), hex2dec('260')]; loop = loop+1;
            
            % Define an Unknown sequence
            tags(loop,:) = [hex2dec('2005'), hex2dec('0010')]; sequences{loop} = [hex2dec('2005'), hex2dec('1085')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1054')]; sequences{loop} = [hex2dec('2005'), hex2dec('1085')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1055')]; sequences{loop} = [hex2dec('2005'), hex2dec('1085')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1056')]; sequences{loop} = [hex2dec('2005'), hex2dec('1085')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1057')]; sequences{loop} = [hex2dec('2005'), hex2dec('1085')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1058')]; sequences{loop} = [hex2dec('2005'), hex2dec('1085')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1059')]; sequences{loop} = [hex2dec('2005'), hex2dec('1085')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('105a')]; sequences{loop} = [hex2dec('2005'), hex2dec('1085')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('105b')]; sequences{loop} = [hex2dec('2005'), hex2dec('1085')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('105c')]; sequences{loop} = [hex2dec('2005'), hex2dec('1085')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('105d')]; sequences{loop} = [hex2dec('2005'), hex2dec('1085')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('105e')]; sequences{loop} = [hex2dec('2005'), hex2dec('1085')]; loop = loop+1;
            
            % Define an Unknown sequence            
            tags(loop,:) = [hex2dec('2005'), hex2dec('0014')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2001'), hex2dec('102d')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2001'), hex2dec('1032')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2001'), hex2dec('1033')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2001'), hex2dec('1035')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2001'), hex2dec('1036')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('0014')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1071')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1072')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1073')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1074')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1075')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1076')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1078')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('1079')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('107a')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('107b')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('107e')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('143c')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('143d')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2005'), hex2dec('143e')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f')]; loop = loop+1;
            % Add the following tags to the base sequence as well
            tags(loop,:) = [hex2dec('2005'), hex2dec('1081')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f'); hex2dec('0'), hex2dec('0')]; loop = loop+1;
            tags(loop,:) = [hex2dec('2001'), hex2dec('0010')]; sequences{loop} = [hex2dec('2001'), hex2dec('105f'); hex2dec('0'), hex2dec('0')]; loop = loop+1;
            
            
            % Define an Unknown sequence
            tags(loop,:) = [hex2dec('0008'), hex2dec('002a')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %DT [20150819]                               #   8,0x 1 AcquisitionDateTime
            tags(loop,:) = [hex2dec('0008'), hex2dec('9123')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %UI [1.3.46.670589.11]                       #  16,0x 1 CreatorVersionUID
            tags(loop,:) = [hex2dec('0008'), hex2dec('9205')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [MONOCHROME]                             #  10,0x 1 PixelPresentation
            tags(loop,:) = [hex2dec('0008'), hex2dec('9206')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [VOLUME]                                 #   6,0x 1 VolumetricProperties
            tags(loop,:) = [hex2dec('0008'), hex2dec('9207')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NONE]                                   #   4,0x 1 VolumeBasedCalculationTechnique
            tags(loop,:) = [hex2dec('0008'), hex2dec('9209')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [T1]                                     #   2,0x 1 AcquisitionContrast
            tags(loop,:) = [hex2dec('0018'), hex2dec('9005')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %SH [T1TFE]                                  #   6,0x 1 PulseSequenceName
            tags(loop,:) = [hex2dec('0018'), hex2dec('9008')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [GRADIENT]                               #   8,0x 1 EchoPulseSequence
            tags(loop,:) = [hex2dec('0018'), hex2dec('9009')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NO]                                     #   2,0x 1 InversionRecovery
            tags(loop,:) = [hex2dec('0018'), hex2dec('9011')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NO]                                     #   2,0x 1 MultipleSpinEcho
            tags(loop,:) = [hex2dec('0018'), hex2dec('9012')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NO]                                     #   2,0x 1 MultiPlanarExcitation
            tags(loop,:) = [hex2dec('0018'), hex2dec('9014')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NO]                                     #   2,0x 1 PhaseContrast
            tags(loop,:) = [hex2dec('0018'), hex2dec('9015')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NO]                                     #   2,0x 1 TimeOfFlightContrast
            tags(loop,:) = [hex2dec('0018'), hex2dec('9016')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [RF]                                     #   2,0x 1 Spoiling
            tags(loop,:) = [hex2dec('0018'), hex2dec('9017')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [LONGITUDINAL]                           #  12,0x 1 SteadyStatePulseSequence
            tags(loop,:) = [hex2dec('0018'), hex2dec('9018')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NO]                                     #   2,0x 1 EchoPlanarPulseSequence
            tags(loop,:) = [hex2dec('0018'), hex2dec('9019')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 7.229999999999997e+75                    #   8,0x 1 TagAngleFirstAxis
            tags(loop,:) = [hex2dec('0018'), hex2dec('9020')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NONE]                                   #   4,0x 1 MagnetizationTransfer
            tags(loop,:) = [hex2dec('0018'), hex2dec('9021')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NO]                                     #   2,0x 1 T2Preparation
            tags(loop,:) = [hex2dec('0018'), hex2dec('9022')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NO]                                     #   2,0x 1 BloodSignalNulling
            tags(loop,:) = [hex2dec('0018'), hex2dec('9024')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NO]                                     #   2,0x 1 SaturationRecovery
            tags(loop,:) = [hex2dec('0018'), hex2dec('9025')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NONE]                                   #   4,0x 1 SpectrallySelectedSuppression
            tags(loop,:) = [hex2dec('0018'), hex2dec('9026')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [WATER]                                  #   6,0x 1 SpectrallySelectedExcitation
            tags(loop,:) = [hex2dec('0018'), hex2dec('9027')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NONE]                                   #   4,0x 1 SpatialPresaturation
            tags(loop,:) = [hex2dec('0018'), hex2dec('9028')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NONE]                                   #   4,0x 1 Tagging
            tags(loop,:) = [hex2dec('0018'), hex2dec('9029')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NONE]                                   #   4,0x 1 OversamplingPhase
            tags(loop,:) = [hex2dec('0018'), hex2dec('9030')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 7.229999999999997e+75                    #   8,0x 1 TagSpacingFirstDimension
            tags(loop,:) = [hex2dec('0018'), hex2dec('9032')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [RECTILINEAR]                            #  12,0x 1 GeometryOfKSpaceTraversal
            tags(loop,:) = [hex2dec('0018'), hex2dec('9033')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [PARTIAL]                                #   8,0x 1 SegmentedKSpaceTraversal
            tags(loop,:) = [hex2dec('0018'), hex2dec('9034')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [UNKNOWN]                                #   8,0x 1 RectilinearPhaseEncodeReordering
            tags(loop,:) = [hex2dec('0018'), hex2dec('9035')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 0                                        #   8,0x 1 TagThickness
            tags(loop,:) = [hex2dec('0018'), hex2dec('9036')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [PHASE]                                  #   6,0x 1 PartialFourierDirection
            tags(loop,:) = [hex2dec('0018'), hex2dec('9037')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NONE]                                   #   4,0x 1 CardiacSynchronizationTechnique
            tags(loop,:) = [hex2dec('0018'), hex2dec('9043')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [MULTICOIL]                              #  10,0x 1 ReceiveCoilType
            tags(loop,:) = [hex2dec('0018'), hex2dec('9044')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NO]                                     #   2,0x 1 QuadratureReceiveCoil
            tags(loop,:) = [hex2dec('0018'), hex2dec('9047')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %SH [MULTI ELEMENT]                          #  14,0x 1 MultiCoilElementName
            tags(loop,:) = [hex2dec('0018'), hex2dec('9050')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %LO					                         #   0,0x 0 TransmitCoilManufacturerName
            tags(loop,:) = [hex2dec('0018'), hex2dec('9051')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [SURFACE]                                #   8,0x 1 TransmitCoilType
            tags(loop,:) = [hex2dec('0018'), hex2dec('9053')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 4.6799999999999997\4.6799999999999997    #  16,0x 2 ChemicalShiftReference
            tags(loop,:) = [hex2dec('0018'), hex2dec('9058')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %US 320                                      #   2,0x 1 MRAcquisitionFrequencyEncodingSteps
            tags(loop,:) = [hex2dec('0018'), hex2dec('9059')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NO]                                     #   2,0x 1 Decoupling
            tags(loop,:) = [hex2dec('0018'), hex2dec('9060')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS											 #   0,0x 0 DecoupledNucleus
            tags(loop,:) = [hex2dec('0018'), hex2dec('9062')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS											 #   0,0x 0 DecouplingMethod
            tags(loop,:) = [hex2dec('0018'), hex2dec('9064')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [RIESZ]                                  #   6,0x 1 KSpaceFiltering
            tags(loop,:) = [hex2dec('0018'), hex2dec('9065')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS											 #   0,0x 0 TimeDomainFiltering
            tags(loop,:) = [hex2dec('0018'), hex2dec('9069')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 1                                        #   8,0x 1 ParallelReductionFactorInPlane
            tags(loop,:) = [hex2dec('0018'), hex2dec('9075')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NONE]                                   #   4,0x 1 DiffusionDirectionality
            tags(loop,:) = [hex2dec('0018'), hex2dec('9077')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NO]                                     #   2,0x 1 ParallelAcquisition
            tags(loop,:) = [hex2dec('0018'), hex2dec('9078')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; % CS										 #   0,0x 0 ParallelAcquisitionTechnique
            tags(loop,:) = [hex2dec('0018'), hex2dec('9079')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 0                                        #   8,0x 1 InversionTimes
            tags(loop,:) = [hex2dec('0018'), hex2dec('9080')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %ST [WATER]                                  #   6,0x 1 MetaboliteMapDescription
            tags(loop,:) = [hex2dec('0018'), hex2dec('9081')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [YES]                                    #   4,0x 1 PartialFourier
            tags(loop,:) = [hex2dec('0018'), hex2dec('9085')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS											 #   0,0x 0 CardiacSignalSource
            tags(loop,:) = [hex2dec('0018'), hex2dec('9090')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 0\0\0                                    #  24,0x 3 VelocityEncodingDirection
            tags(loop,:) = [hex2dec('0018'), hex2dec('9091')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 0                                        #   8,0x 1 VelocityEncodingMinimumValue
            tags(loop,:) = [hex2dec('0018'), hex2dec('9093')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %US 1                                        #   2,0x 1 NumberOfKSpaceTrajectories
            tags(loop,:) = [hex2dec('0018'), hex2dec('9094')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS											 #   0,0x 0 CoverageOfKSpace
            tags(loop,:) = [hex2dec('0018'), hex2dec('9101')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NO]                                     #   2,0x 1 FrequencyCorrection
            tags(loop,:) = [hex2dec('0018'), hex2dec('9147')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS											 #   0,0x 0 DiffusionAnisotropyType
            tags(loop,:) = [hex2dec('0018'), hex2dec('9155')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 1                                        #   8,0x 1 ParallelReductionFactorOutOfPlane
            tags(loop,:) = [hex2dec('0018'), hex2dec('9168')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 1                                        #   8,0x 1 ParallelReductionFactorSecondInPlane
            tags(loop,:) = [hex2dec('0018'), hex2dec('9169')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS											 #   0,0x 0 CardiacBeatRejectionTechnique
            tags(loop,:) = [hex2dec('0018'), hex2dec('9170')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NONE]                                   #   4,0x 1 RespiratoryMotionCompensationTechnique
            tags(loop,:) = [hex2dec('0018'), hex2dec('9171')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NONE]                                   #   4,0x 1 RespiratorySignalSource
            tags(loop,:) = [hex2dec('0018'), hex2dec('9172')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NONE]                                   #   4,0x 1 BulkMotionCompensationTechnique
            tags(loop,:) = [hex2dec('0018'), hex2dec('9174')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [IEC]                                    #   4,0x 1 ApplicableSafetyStandardAgency
            tags(loop,:) = [hex2dec('0018'), hex2dec('9177')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [STATIC FIELD]                           #  12,0x 1 OperatingModeType
            tags(loop,:) = [hex2dec('0018'), hex2dec('9178')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [IEC_FIRST_LEVEL]                        #  16,0x 1 OperatingMode
            tags(loop,:) = [hex2dec('0018'), hex2dec('9177')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [RF]                                     #   2,0x 1 OperatingModeType
            tags(loop,:) = [hex2dec('0018'), hex2dec('9178')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [IEC_NORMAL]                             #  10,0x 1 OperatingMode
            tags(loop,:) = [hex2dec('0018'), hex2dec('9177')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [GRADIENT]                               #   8,0x 1 OperatingModeType
            tags(loop,:) = [hex2dec('0018'), hex2dec('9178')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [IEC_NORMAL]                             #  10,0x 1 OperatingMode
            tags(loop,:) = [hex2dec('0018'), hex2dec('9179')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [IEC_WHOLE_BODY]                         #  14,0x 1 SpecificAbsorptionRateDefinition
            tags(loop,:) = [hex2dec('0018'), hex2dec('9180')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [DB_DT]                                  #   6,0x 1 GradientOutputType
            tags(loop,:) = [hex2dec('0018'), hex2dec('9181')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 0.38351830840110779                      #   8,0x 1 SpecificAbsorptionRateValue
            tags(loop,:) = [hex2dec('0018'), hex2dec('9182')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 35.410717010498047                       #   8,0x 1 GradientOutput
            tags(loop,:) = [hex2dec('0018'), hex2dec('9183')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS											 #   0,0x 0 FlowCompensationDirection
            tags(loop,:) = [hex2dec('0018'), hex2dec('9199')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [NO]                                     #   2,0x 1 WaterReferencedPhaseCorrection
            tags(loop,:) = [hex2dec('0018'), hex2dec('9200')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS											 #   0,0x 0 MRSpectroscopyAcquisitionType
            tags(loop,:) = [hex2dec('0018'), hex2dec('9231')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %US 214                                      #   2,0x 1 MRAcquisitionPhaseEncodingStepsInPlane
            tags(loop,:) = [hex2dec('0018'), hex2dec('9232')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %US 1                                        #   2,0x 1 MRAcquisitionPhaseEncodingStepsOutOfPlane
            tags(loop,:) = [hex2dec('0018'), hex2dec('9240')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %US 0                                        #   2,0x 1 RFEchoTrainLength
            tags(loop,:) = [hex2dec('0018'), hex2dec('9241')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %US 97                                       #   2,0x 1 GradientEchoTrainLength
            tags(loop,:) = [hex2dec('0018'), hex2dec('9602')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 7.229999999999997e+75                    #   8,0x 1 DiffusionBValueXX
            tags(loop,:) = [hex2dec('0018'), hex2dec('9603')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 7.229999999999997e+75                    #   8,0x 1 DiffusionBValueXY
            tags(loop,:) = [hex2dec('0018'), hex2dec('9604')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 7.229999999999997e+75                    #   8,0x 1 DiffusionBValueXZ
            tags(loop,:) = [hex2dec('0018'), hex2dec('9605')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 7.229999999999997e+75                    #   8,0x 1 DiffusionBValueYY
            tags(loop,:) = [hex2dec('0018'), hex2dec('9606')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 7.229999999999997e+75                    #   8,0x 1 DiffusionBValueYZ
            tags(loop,:) = [hex2dec('0018'), hex2dec('9607')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 7.229999999999997e+75                    #   8,0x 1 DiffusionBValueZZ
            tags(loop,:) = [hex2dec('0020'), hex2dec('9072')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS [U]                                      #   2,0x 1 FrameLaterality
            tags(loop,:) = [hex2dec('0020'), hex2dec('9254')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 0                                        #   8,0x 1 RespiratoryIntervalTime
            tags(loop,:) = [hex2dec('0020'), hex2dec('9255')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 0                                        #   8,0x 1 NominalRespiratoryTriggerDelayTime
            tags(loop,:) = [hex2dec('0020'), hex2dec('9256')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %FD 7.229999999999997e+75                    #   8,0x 1 RespiratoryTriggerDelayThreshold
            tags(loop,:) = [hex2dec('0028'), hex2dec('9001')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %UL 1                                        #   4,0x 1 DataPointRows
            tags(loop,:) = [hex2dec('0028'), hex2dec('9002')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %UL 0                                        #   4,0x 1 DataPointColumns
            tags(loop,:) = [hex2dec('0028'), hex2dec('9003')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS											 #   0,0x 0 SignalDomainColumns
            tags(loop,:) = [hex2dec('0028'), hex2dec('9108')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %CS											 #   0,0x 0 DataRepresentation
            tags(loop,:) = [hex2dec('0040'), hex2dec('9210')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f')]; loop = loop+1; %SH [Philips]                                #   8,0x 1 LUTLabel
            
            % Add the following tags to the base sequence as well
            tags(loop,:) = [hex2dec('0018'),hex2dec('9087')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f'); hex2dec('0'), hex2dec('0')]; loop = loop+1; %FD 7.229999999999997e+75                    #   8,0x 1 DiffusionBValue
            tags(loop,:) = [hex2dec('0018'),hex2dec('9089')]; sequences{loop} = [hex2dec('2005'),hex2dec('140f'); hex2dec('0'), hex2dec('0')]; loop = loop+1; %FD 7.229999999999997e+75\7.229999999999997e+75\7.229999999999997e+75 #  24,0x 3 DiffusionGradientOrientation
            
            D.mSequences = {tags, sequences};
        end
    end
    methods (Static, Hidden, Access = protected )
        function date_time_vec = format_date_time(date_time)
            date_time_vec = {};
            splitted_string = strsplit(date_time, '/');
            if( length(splitted_string) ~= 2 )
                return;
            end
            splitted_date = strsplit(splitted_string{1}, '.');
            if( length(splitted_date) == 3 )
                date_time_vec{1} = [strtrim(splitted_date{3}), strtrim(splitted_date{2}), strtrim(splitted_date{1})];
            end
            splitted_time = strsplit(splitted_string{2}, ':');
            if( length(splitted_time) >= 3 )
                date_time_vec{2} = [strtrim(splitted_time{1}), strtrim(splitted_time{2}), strtrim(splitted_time{3}), '.00000'];
            end
        end
        function uid = get_uid( base_uid )
            temp_uid = dicomuid;
            
            % find the last dot
            dotind = strfind( base_uid, '.' );
            
            dotind_temp = strfind( temp_uid, '.');
            if ~isempty(dotind) && ~isempty(temp_uid)
                dotind = dotind(end);
                dotind_temp = dotind_temp(end);
                uid = base_uid(1:dotind);
                uid = [uid, temp_uid(dotind_temp+1:end)];
            else
                uid = dicomuid;
            end
            
            if( length(uid) > length(base_uid) )
                uid = uid(1:length(base_uid));
            elseif( length(uid) < length(base_uid) )
                uid(end+1:length(base_uid)) = '0';
            end
        end
        function age = get_patients_age(patients_birth_date)
            current_year = str2double(datestr(now,'yyyy'));
            
            % extract year from birth date
            if (length(patients_birth_date) > 3)
                year = str2double(patients_birth_date(1:4));
                age = [num2str(current_year - year), 'Y'];
            end
        end
        function [aDirectionCosinesRL, aDirectionCosinesAP, aDirectionCosinesFH] = get_direction_cosines(Parfile, image_nr)
            vRightHanded = true;
            
            vDegToRad = 0.01745329252;
            
            vAngle = zeros(3, 1);
            
            switch (Parfile.ImageInformation.SliceOrientation(image_nr))
                case 3 %COR:
                    vViewMatrixX = zeros(3, 1); vViewMatrixX(1) = 1.0;
                    vViewMatrixY = zeros(3, 1); vViewMatrixY(3) = 1.0;
                    vViewMatrixZ = zeros(3, 1); vViewMatrixZ(2) = -1.0;
                case 2 %SAG:
                    vViewMatrixX = zeros(3, 1); vViewMatrixX(3) = -1.0;
                    vViewMatrixY = zeros(3, 1); vViewMatrixY(1) = 1.0;
                    vViewMatrixZ = zeros(3, 1); vViewMatrixZ(2) = -1.0;
                case 1 %TRA:
                    vViewMatrixX = zeros(3, 1); vViewMatrixX(1) = 1.0;
                    vViewMatrixY = zeros(3, 1); vViewMatrixY(2) = 1.0;
                    vViewMatrixZ = zeros(3, 1); vViewMatrixZ(3) = 1.0;
                otherwise
                    vViewMatrixX = zeros(3, 1); vViewMatrixX(1) = 1.0;
                    vViewMatrixY = zeros(3, 1); vViewMatrixY(2) = 1.0;
                    vViewMatrixZ = zeros(3, 1); vViewMatrixZ(3) = 1.0;
            end
            
            vAngle(1) = Parfile.ImageInformation.ImageAngulation(image_nr, 3);
            vAngle(2) = Parfile.ImageInformation.ImageAngulation(image_nr, 1);
            vAngle(3) = Parfile.ImageInformation.ImageAngulation(image_nr, 2);
            
            
            sx = sin(vAngle(1) * vDegToRad);
            sy = sin(vAngle(2) * vDegToRad);
            sz = sin(vAngle(3) * vDegToRad);
            cx = cos(vAngle(1) * vDegToRad);
            cy = cos(vAngle(2) * vDegToRad);
            cz = cos(vAngle(3) * vDegToRad);
            
            if (vRightHanded)
                sx = -sx;
                sy = -sy;
                sz = -sz;
            end
            
            vRotMatrixX = zeros(3, 1);
            vRotMatrixY = zeros(3, 1);
            vRotMatrixZ = zeros(3, 1);
            
            vRotMatrixX(1) = cy * cz;
            vRotMatrixY(1) = -sz * cx + sx * sy * cz;
            vRotMatrixZ(1) = sx * sz + sy * cx * cz;
            
            vRotMatrixX(2) = sz * cy;
            vRotMatrixY(2) = cx * cz + sx * sy * sz;
            vRotMatrixZ(2) = -sx * cz + sy * sz * cx;
            
            vRotMatrixX(3) = -sy;
            vRotMatrixY(3) = sx * cy;
            vRotMatrixZ(3) = cx * cy;
            
            aDirectionCosinesRL = zeros(3, 1);
            aDirectionCosinesAP = zeros(3, 1);
            aDirectionCosinesFH = zeros(3, 1);
            
            aDirectionCosinesRL(1) = vRotMatrixX(1) * vViewMatrixX(1) + vRotMatrixX(2) * vViewMatrixY(1) + vRotMatrixX(3) * vViewMatrixZ(1);
            aDirectionCosinesRL(2) = vRotMatrixX(1) * vViewMatrixX(2) + vRotMatrixX(2) * vViewMatrixY(2) + vRotMatrixX(3) * vViewMatrixZ(2);
            aDirectionCosinesRL(3) = vRotMatrixX(1) * vViewMatrixX(3) + vRotMatrixX(2) * vViewMatrixY(3) + vRotMatrixX(3) * vViewMatrixZ(3);
            aDirectionCosinesAP(1) = vRotMatrixY(1) * vViewMatrixX(1) + vRotMatrixY(2) * vViewMatrixY(1) + vRotMatrixY(3) * vViewMatrixZ(1);
            aDirectionCosinesAP(2) = vRotMatrixY(1) * vViewMatrixX(2) + vRotMatrixY(2) * vViewMatrixY(2) + vRotMatrixY(3) * vViewMatrixZ(2);
            aDirectionCosinesAP(3) = vRotMatrixY(1) * vViewMatrixX(3) + vRotMatrixY(2) * vViewMatrixY(3) + vRotMatrixY(3) * vViewMatrixZ(3);
            aDirectionCosinesFH(1) = vRotMatrixZ(1) * vViewMatrixX(1) + vRotMatrixZ(2) * vViewMatrixY(1) + vRotMatrixZ(3) * vViewMatrixZ(1);
            aDirectionCosinesFH(2) = vRotMatrixZ(1) * vViewMatrixX(2) + vRotMatrixZ(2) * vViewMatrixY(2) + vRotMatrixZ(3) * vViewMatrixZ(2);
            aDirectionCosinesFH(3) = vRotMatrixZ(1) * vViewMatrixX(3) + vRotMatrixZ(2) * vViewMatrixY(3) + vRotMatrixZ(3) * vViewMatrixZ(3);
        end
        function item = sort_tags(item)
            fn = fieldnames(item);
            GroupElement = zeros(length(fn),2);
            
            for i = 1:length(fn)
                [g,e] = dicomlookup(fn{i});
                if( isempty(g) )
                    GroupElement(i,1) = 65535;
                else
                    GroupElement(i,1) = g;
                end
                if( isempty(e) )
                    GroupElement(i,2) = 65535;
                else
                    GroupElement(i,2) = e;
                end
                
                if( isstruct( item.(fn{i}) ) )
                    fn2 = fieldnames(item.(fn{i}));
                    if( strcmpi(fn2{1}, 'Item_1' ))
                        item.(fn{i}).Item_1 = DICOMExporter.sort_tags(item.(fn{i}).Item_1);
                    end
                end
            end
            [temp,ind] = sortrows(GroupElement);
            item = orderfields(item, ind);
        end
        function create_dicomdir( out_dir, files_dir )
            % use the external too l from dcmtk to create the DICOMDIR file
            % (if available)
            base_path = which('MRecon.m');
            if ~isempty( base_path )
                base_path = strrep( base_path, 'MRecon.m', '' );
            else
                base_path = which('MRecon.p');
                if ~isempty( base_path )
                    base_path = strrep( base_path, 'MRecon.p', '' );
                end
            end
            if ~isempty(base_path) && exist( [base_path,'/tools'], 'dir' )
                tools_path = [base_path, filesep, 'tools'];
                dicomdir_file = [out_dir, filesep, 'DICOMDIR'];
                if( exist(dicomdir_file, 'file' ))
                    delete( dicomdir_file );
                end
                try
                    arguments = ['+D ', out_dir, filesep, 'DICOMDIR +id ', files_dir, ' +r'];
                    if ispc
                        system( [tools_path, filesep, 'dcmmkdir.exe ', arguments] );
                    elseif ismac
                        system( [tools_path, filesep, 'dcmmkdir_mac ', arguments] );
                    elseif isunix
                        system( [tools_path, filesep, 'dcmmkdir_linux ', arguments] );
                    end
                catch
                    warning('Could not create the DICOMDIR file');
                end
            end
        end      
    end
end
