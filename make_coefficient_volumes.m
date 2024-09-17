function make_coefficient_volumes()
    % by Saskia. Takes collated results from visualisation jobs and creates
    % volumes, ready to be projected onto the cortical surface.
    
    addpath('/group/mlr-lab/AH/Projects/spm12/')
    addpath('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/')
    root = '/imaging/projects/cbu/wbic-p00591-DAISY/main/work/';
    cd(root);

    % load template volume - in MNI space, EPI resolution
    vol = spm_vol('/imaging/projects/cbu/wbic-p00591-DAISY/main/work/sub-001/run-01/tmpT1.nii.gz');
    T = vol.mat;
    template = spm_read_vols(vol);
    
    subcode = {'sub-001','sub-002','sub-003','sub-004','sub-007','sub-009','sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019','sub-020','sub-021','sub-022','sub-023','sub-024','sub-026','sub-028','sub-029','sub-030','sub-031','sub-032'};
    
    % get metadata
    load('/imaging/projects/cbu/wbic-p00591-DAISY/main/derivatives/cox/tedana_averaged_metadata.mat')

    % get each participant's coordinates and put them into MNI space
    mni = {};
    for s = 1:size(subcode,2)
        mni{s,1} = mni2cor(metadata(s).coords.xyz,T);
    end
    
    % 1. CLASSIFICATION (Logistic LASSO and SOSLASSO)
    
    % 1A. Logistic LASSO
    
    % for every participant
    for s = 1:size(subcode,2)
    
        % make output directory
        delete([root,'/',subcode{s},'/coefficients/surface/*'])
        delete([root,'/',subcode{s},'/coefficients/volume/*'])
        mkdir([root,'/',subcode{s},'/coefficients/volume/'])

        % load results   
        load([root,'/',subcode{s},'/coefficients/mat/tedana/log.mat']);

        % get coordinates
        coords = mni{s,1};

        % initialise volume
        volume = zeros(size(template));

        % set voxels to final coefficient values
        for v = 1:size(coords,1)
            volume(coords(v,1),coords(v,2),coords(v,3)) = finalcoefs(v,1);
        end

        % save 
        vol.fname = [root,'/',subcode{s},'/coefficients/volume/',subcode{s},'_log-LASSO_final_coefficients.nii'];
        spm_write_vol(vol,volume);
        
        % for every perm, do the same
        for perm = 1:100        
            volume = zeros(size(template));
            for v = 1:size(coords,1)
                volume(coords(v,1),coords(v,2),coords(v,3)) = permcoefs(v,perm);
            end
            vol.fname = [root,'/',subcode{s},'/coefficients/volume/',subcode{s},'_log-LASSO_randomseed-',num2str(perm),'_perm_coefficients.nii'];
            spm_write_vol(vol,volume);           
        end 
    end
    
    % 1B. SOSLASSO
    
    % load results from condor
    
    % load results
    load(['/imaging/projects/cbu/wbic-p00591-DAISY/main/derivatives/condor/derivatives/classification/SOSLASSO/visualize/final_performance.mat']);
    
    % for each participant
    for s = 1:size(subcode,2)

        % get just the row applying to that participant
        tmp = str2num(erase(subcode{s},'sub-'));
        subtable = Tallcv(Tallcv.subject==tmp,:);

        % get coordinates
        coords = mni{s,1};

        % initialise volume
        volume = zeros(size(template));

        % set voxels to final coefficient values
        finalcoefs = subtable{1,1}{:,:};
        for v = 1:size(coords,1)
            volume(coords(v,1),coords(v,2),coords(v,3)) = finalcoefs(v,1);
        end

        % save 
        vol.fname = [root,'/',subcode{s},'/coefficients/volume/',subcode{s},'_SOSLASSO_final_coefficients.nii'];
        spm_write_vol(vol,volume);  
    end

    % load perm results
    load(['/imaging/projects/cbu/wbic-p00591-DAISY/main/derivatives/condor/derivatives/classification/SOSLASSO/visualize/perm_performance.mat']);
    
    % for each participant
    for s = 1:size(subcode,2)

        % get just the rows applying to that participant
        tmp = str2num(erase(subcode{s},'sub-'));
        subtable = Tallcv(Tallcv.subject==tmp,:);

        % get coordinates
        coords = mni{s,1};
        
        % for every random seed
        for perm = 1:100

            subsubtable = subtable(subtable.RandomSeed == perm,:);
            
            % initialise volume
            volume = zeros(size(template));

            % set voxels to final coefficient values
            permcoefs = subsubtable{1,1}{:,:};
            for v = 1:size(coords,1)
                volume(coords(v,1),coords(v,2),coords(v,3)) = permcoefs(v,1);
            end

            % save
            vol.fname = [root,'/',subcode{s},'/coefficients/volume/',subcode{s},'_SOSLASSO_randomseed-',num2str(perm),'_perm_coefficients.nii'];
            spm_write_vol(vol,volume);
        end
    end

    % 2. CORRELATION (RSL with linear LASSO and with grOWL)

    
    % 2A. Linear LASSO

    % for every participant
    for s = 1:size(subcode,2)
    
        % for every dimension
        for d = 1:3

            % load results   
            load([root,'/',subcode{s},'/coefficients/mat/tedana/dimension_',num2str(d),'_lin.mat']);
    
            % get coordinates
            coords = mni{s,1};
    
            % initialise volume
            volume = zeros(size(template));
    
            % set voxels to final coefficient values
            for v = 1:size(coords,1)
                volume(coords(v,1),coords(v,2),coords(v,3)) = finalcoefs(v,1);
            end
    
            % save 
            vol.fname = [root,'/',subcode{s},'/coefficients/volume/',subcode{s},'_linear-LASSO_dimension-',num2str(d),'_final_coefficients.nii'];
            spm_write_vol(vol,volume);
            
            % for every perm, do the same
            for perm = 1:100        
                volume = zeros(size(template));
                for v = 1:size(coords,1)
                    volume(coords(v,1),coords(v,2),coords(v,3)) = permcoefs(v,perm);
                end
                vol.fname = [root,'/',subcode{s},'/coefficients/volume/',subcode{s},'_linear-LASSO_dimension-',num2str(d),'_randomseed-',num2str(perm),'_perm_coefficients.nii'];
                spm_write_vol(vol,volume);           
            end 
        end
    end
    
    % 1B. grOWL
    
    % load results from condor
    
    % load results
    load(['/imaging/projects/cbu/wbic-p00591-DAISY/main/derivatives/condor/derivatives/correlation/grOWL/visualize/wholebrain/final_performance.mat']);
    
    % for each participant
    for s = 1:size(subcode,2)

        % get just the row applying to that participant
        tmp = str2num(erase(subcode{s},'sub-'));
        subtable = Tallcv(Tallcv.subject==tmp,:);

        % get coordinates. These coefficients are sparse so we must index
        % the coefficients that we need
        coords = mni{s,1};
        coords = coords(subtable{1,21}{:,:},:);

        % set voxels to final coefficient values
        finalcoefs = subtable{1,1}{:,:};

        % for each dimension
        for d = 1:3
            
            % initialise volume
            volume = zeros(size(template));

            for v = 1:size(coords,1)
                volume(coords(v,1),coords(v,2),coords(v,3)) = finalcoefs(v,d);
            end
    
            % save 
            vol.fname = [root,'/',subcode{s},'/coefficients/volume/',subcode{s},'_grOWL_dimension-',num2str(d),'_final_coefficients.nii'];
            spm_write_vol(vol,volume);  
        end
    end

    % load perm results
    load(['/imaging/projects/cbu/wbic-p00591-DAISY/main/derivatives/condor/derivatives/correlation/grOWL/visualize/wholebrain/perm_performance.mat']);
    
    % for each participant
    for s = 1:size(subcode,2)

        % get just the rows applying to that participant
        tmp = str2num(erase(subcode{s},'sub-'));
        subtable = Tallcv(Tallcv.subject==tmp,:);
        
        % for every random seed
        for perm = 1:100

            subsubtable = subtable(subtable.RandomSeed == perm,:);
            
            % for each dimension 
            for d = 1:3

                % get coordinates. Again, these are sparse
                coords = mni{s,1};
                coords = coords(subsubtable{1,21}{:,:},:);

                % initialise volume
                volume = zeros(size(template));
    
                % set voxels to final coefficient values
                permcoefs = subsubtable{1,1}{:,:};
                for v = 1:size(coords,1)
                    volume(coords(v,1),coords(v,2),coords(v,3)) = permcoefs(v,d);
                end
    
                % save
                vol.fname = [root,'/',subcode{s},'/coefficients/volume/',subcode{s},'_grOWL_dimension-',num2str(d),'_randomseed-',num2str(perm),'_perm_coefficients.nii'];
                spm_write_vol(vol,volume);
            end
        end
    end
    