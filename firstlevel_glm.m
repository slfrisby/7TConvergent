function firstlevel_glm(ids)

    addpath('/group/mlr-lab/AH/Projects/spm12/');
    addpath('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts');

    % % if design matrices do not exist, construct them. Do not do this
    % on the submit node
    % if ~exist('/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices')
    %     make_design_matrices
    % end
    % clear subcode

    root = ['/imaging/projects/cbu/wbic-p00591-DAISY/main/derivatives'];
    cd([root]);
    
    subcode{1} = ['sub-',ids] 

    runs = {'01','02','03','04'};

    % for every participant
    for s = 1:size(subcode,2)

        cd(root)
        % check that motion files exist for all 4 runs
        if size(dir([root,'/halaiprep/',subcode{s},'/func/',subcode{s},'_*_motion.txt']),1) == 4
            disp('*** FILES FOUND - PROCESSING ***')
            % make output directories 
            mkdir([root,'/GLM/first/native/t2star/',subcode{s}])
            mkdir([root,'/GLM/first/native/tedana/',subcode{s}])
            mkdir([root,'/GLM/first/mni/t2star/',subcode{s}])
            mkdir([root,'/GLM/first/mni/tedana/',subcode{s}])
            mkdir([root,'/SPM/first/native/',subcode{s}])
            mkdir([root,'/SPM/first/mni/',subcode{s}])
            % if SPM directories already exist, empty them
            delete([root,'/SPM/first/native/',subcode{s},'/*'])
            delete([root,'/SPM/first/mni/',subcode{s},'/*'])

            % unzip native-space grey matter mask into SPM directory
            gunzip([root,'/fmriprep/',subcode{s},'/anat/',subcode{s},'_run-1_space-native_label-GM_mask.nii.gz'],[root,'/SPM/first/native/',subcode{s},'/']);

            % for each run
            for r = 1:length(runs)

                % load motion file and extract 6 motion parameters. These
                % are the same for both native and MNI space models.
                R = spm_load([root,'/halaiprep/',subcode{s},'/func/',subcode{s},'_run-',runs{r},'_motion.txt']);
                save([root,'/SPM/first/native/',subcode{s},'/motion_run-',runs{r},'.mat'],'R'); 
                save([root,'/SPM/first/mni/',subcode{s},'/motion_run-',runs{r},'.mat'],'R'); 

                % unzip t2star functional files
                tmp=dir(['./halaiprep/',subcode{s},'/func/',subcode{s},'_*_rec-t2star_run-',runs{r},'_space-native_desc-preproc_bold.nii.gz']);
                [~,name,~]=fileparts(tmp(1).name);
                gunzip([tmp(1).folder,'/',tmp(1).name],[root,'/SPM/first/native/',subcode{s},'/']);
                tmp=dir(['./halaiprep/',subcode{s},'/func/',subcode{s},'_*_rec-t2star_run-',runs{r},'_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz']);
                [~,name,~]=fileparts(tmp(1).name);
                gunzip([tmp(1).folder,'/',tmp(1).name],[root,'/SPM/first/mni/',subcode{s},'/']);

                % create smoothed version for use in univariate analysis.
                clear matlabbatch
                matlabbatch{1}.spm.spatial.smooth.data = cellstr(spm_select('ExtFPList',[root,'/SPM/first/mni/',subcode{s},'/'],name,[1:551]));
                matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];% in mm. This is pre-specified by when the function is called
                matlabbatch{1}.spm.spatial.smooth.dtype = 0; % specifies what sort of output file to write (0 = same as input)
                matlabbatch{1}.spm.spatial.smooth.im = 0; % don't use implicit mask (which assumes that any voxels set to 0 or NaN in the input should stay the same in the output)
                matlabbatch{1}.spm.spatial.smooth.prefix = ['s6'];
                spm_jobman('run',matlabbatch);

                % repeat for tedana. Unzip functional files
                tmp=dir(['./halaiprep/',subcode{s},'/func/',subcode{s},'_*_rec-tedana_run-',runs{r},'_space-native_desc-preproc_bold.nii.gz']);
                [~,name,~]=fileparts(tmp(1).name);
                gunzip([tmp(1).folder,'/',tmp(1).name],[root,'/SPM/first/native/',subcode{s},'/']);
                tmp=dir(['./halaiprep/',subcode{s},'/func/',subcode{s},'_*_rec-tedana_run-',runs{r},'_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz']);
                [~,name,~]=fileparts(tmp(1).name);
                gunzip([tmp(1).folder,'/',tmp(1).name],[root,'/SPM/first/mni/',subcode{s},'/']);

                % create smoothed version
                clear matlabbatch
                matlabbatch{1}.spm.spatial.smooth.data = cellstr(spm_select('ExtFPList',[root,'/SPM/first/mni/',subcode{s},'/'],name,[1:551]));
                matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];% in mm. This is pre-specified by when the function is called
                matlabbatch{1}.spm.spatial.smooth.dtype = 0; % specifies what sort of output file to write (0 = same as input)
                matlabbatch{1}.spm.spatial.smooth.im = 0; % don't use implicit mask (which assumes that any voxels set to 0 or NaN in the input should stay the same in the output)
                matlabbatch{1}.spm.spatial.smooth.prefix = ['s6'];
                spm_jobman('run',matlabbatch);

            end

            %% fit first-level GLM for decoding (native space)
            
            % t2star - build GLM
            clear matlabbatch
            matlabbatch{1}.spm.stats.fmri_spec.dir = {[root,'/GLM/first/native/t2star/',subcode{s},'/']}; % directory for outputting full design matrix (this means the final one for use in the GLM, i.e. the bit specified manually plus motion etc.)
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % units in which to provide timing for the parts of the design matrix that are specified manually
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.51; % TR
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 48; % number of slices
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 48/2; % align predictor variables so they are predicting responses midway through volume acquisition
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/native/',subcode{s},'/'],[subcode{s},'_task-namingslowish_acq-MEMB_rec-t2star_run-01_space-native_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {[root,'/SPM/first/native/',subcode{s},'/motion_run-01.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-01_type-multivariate_design-matrix.mat']}; % design matrix specific to this run
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128; % high-pass filter - slow signal drifts with a period LONGER than this will be removed (this is the default value)
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {}); % format to expect design matrix in
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {}); % format to expect motion regressors in
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/native/',subcode{s},'/'],[subcode{s},'_task-namingslowish_acq-MEMB_rec-t2star_run-02_space-native_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi_reg = {[root,'/SPM/first/native/',subcode{s},'/motion_run-02.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-02_type-multivariate_design-matrix.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).hpf = 128;
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/native/',subcode{s},'/'],[subcode{s},'_task-namingslowish_acq-MEMB_rec-t2star_run-03_space-native_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).multi_reg = {[root,'/SPM/first/native/',subcode{s},'/motion_run-03.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-03_type-multivariate_design-matrix.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).hpf = 128; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/native/',subcode{s},'/'],[subcode{s},'_task-namingslowish_acq-MEMB_rec-t2star_run-04_space-native_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).multi_reg = {[root,'/SPM/first/native/',subcode{s},'/motion_run-04.mat']};
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-04_type-multivariate_design-matrix.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).hpf = 128; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.mask = {[root,'/SPM/first/native/',subcode{s},'/',subcode{s},'_run-1_space-native_label-GM_mask.nii']}; % native space grey matter mask
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; % do not assume that the HRF varies between subjects and voxels
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1; % use linear convolution for the HRF (https://pubmed.ncbi.nlm.nih.gov/9438436/)
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None'; % no global intensity normalisation
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0; % default to go with the above
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; % model autocorrelation

            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; % as opposed to Bayesian
            spm_jobman('run',matlabbatch);

            % tedana - build GLM
            clear matlabbatch
            matlabbatch{1}.spm.stats.fmri_spec.dir = {[root,'/GLM/first/native/tedana/',subcode{s},'/']}; % directory for outputting full design matrix (this means the final one for use in the GLM, i.e. the bit specified manually plus motion etc.)
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % units in which to provide timing for the parts of the design matrix that are specified manually
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.51; % TR
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 48; % number of slices
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 48/2; % align predictor variables so they are predicting responses midway through volume acquisition
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/native/',subcode{s},'/'],[subcode{s},'_task-namingslowish_acq-MEMB_rec-tedana_run-01_space-native_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {[root,'/SPM/first/native/',subcode{s},'/motion_run-01.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-01_type-multivariate_design-matrix.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128; % high-pass filter - slow signal drifts with a period LONGER than this will be removed (this is the default value)
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {}); % format to expect design matrix in
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {}); % format to expect motion regressors in
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/native/',subcode{s},'/'],[subcode{s},'_task-namingslowish_acq-MEMB_rec-tedana_run-02_space-native_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi_reg = {[root,'/SPM/first/native/',subcode{s},'/motion_run-02.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-02_type-multivariate_design-matrix.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).hpf = 128;
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/native/',subcode{s},'/'],[subcode{s},'_task-namingslowish_acq-MEMB_rec-tedana_run-03_space-native_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).multi_reg = {[root,'/SPM/first/native/',subcode{s},'/motion_run-03.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-03_type-multivariate_design-matrix.mat']};
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).hpf = 128; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/native/',subcode{s},'/'],[subcode{s},'_task-namingslowish_acq-MEMB_rec-tedana_run-04_space-native_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).multi_reg = {[root,'/SPM/first/native/',subcode{s},'/motion_run-04.mat']};
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-04_type-multivariate_design-matrix.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).hpf = 128; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.mask = {[root,'/SPM/first/native/',subcode{s},'/',subcode{s},'_run-1_space-native_label-GM_mask.nii']}; % native space grey matter mask
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; % do not assume that the HRF varies between subjects and voxels
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1; % use linear convolution for the HRF (https://pubmed.ncbi.nlm.nih.gov/9438436/)
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None'; % no global intensity normalisation
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0; % default to go with the above
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; % model autocorrelation

            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; % as opposed to Bayesian
            spm_jobman('run',matlabbatch);

            %% fit first-level GLM for univariate analysis (MNI space)

            % t2star - build GLM
            clear matlabbatch
            matlabbatch{1}.spm.stats.fmri_spec.dir = {[root,'/GLM/first/mni/t2star/',subcode{s},'/']}; % directory for outputting full design matrix (this means the final one for use in the GLM, i.e. the bit specified manually plus motion etc.)
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % units in which to provide timing for the parts of the design matrix that are specified manually
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.51; % TR
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 48; % number of slices
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 48/2; % align predictor variables so they are predicting responses midway through volume acquisition
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/mni/',subcode{s},'/'],['s6',subcode{s},'_task-namingslowish_acq-MEMB_rec-t2star_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {[root,'/SPM/first/mni/',subcode{s},'/motion_run-01.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-01_type-univariate_design-matrix.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128; % high-pass filter - slow signal drifts with a period LONGER than this will be removed (this is the default value)
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {}); % format to expect design matrix in
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {}); % format to expect motion regressors in
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/mni/',subcode{s},'/'],['s6',subcode{s},'_task-namingslowish_acq-MEMB_rec-t2star_run-02_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi_reg = {[root,'/SPM/first/mni/',subcode{s},'/motion_run-02.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-02_type-univariate_design-matrix.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).hpf = 128;
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/mni/',subcode{s},'/'],['s6',subcode{s},'_task-namingslowish_acq-MEMB_rec-t2star_run-03_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).multi_reg = {[root,'/SPM/first/mni/',subcode{s},'/motion_run-03.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-03_type-univariate_design-matrix.mat']};
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).hpf = 128; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/mni/',subcode{s},'/'],['s6',subcode{s},'_task-namingslowish_acq-MEMB_rec-t2star_run-04_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).multi_reg = {[root,'/SPM/first/mni/',subcode{s},'/motion_run-04.mat']};
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-04_type-univariate_design-matrix.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).hpf = 128; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.mask = {'/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'}; % MNI space brain mask
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; % do not assume that the HRF varies between subjects and voxels
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1; % use linear convolution for the HRF (https://pubmed.ncbi.nlm.nih.gov/9438436/)
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None'; % no global intensity normalisation
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0; % default to go with the above
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; % model autocorrelation

            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; % as opposed to Bayesian

            matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'Task>Fix';
            matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = [[1 -1; -1 1] zeros(2,8) [1 -1; -1 1] zeros(2,8) [1 -1; -1 1] zeros(2,8) [1 -1; -1 1] zeros(2,8)];
            matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'A';
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 zeros(1,9) 1 zeros(1,9) 1 zeros(1,9) 1 zeros(1,9)];
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'I';
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 1 zeros(1,9) 1 zeros(1,9) 1 zeros(1,9) 1 zeros(1,8)];
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'A>I';
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [1 -1 zeros(1,8) 1 -1 zeros(1,8) 1 -1 zeros(1,8) 1 -1 zeros(1,8)];
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'I>A';
            matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [-1 1 zeros(1,8) -1 1 zeros(1,8) -1 1 zeros(1,8) -1 1 zeros(1,8) ];
            matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'A+I';
            matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [1 1 zeros(1,8) 1 1 zeros(1,8) 1 1 zeros(1,8) 1 1 zeros(1,8) ];
            matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
            spm_jobman('run',matlabbatch);

            % tedana - build GLM
            clear matlabbatch
            matlabbatch{1}.spm.stats.fmri_spec.dir = {[root,'/GLM/first/mni/tedana/',subcode{s},'/']}; % directory for outputting full design matrix (this means the final one for use in the GLM, i.e. the bit specified manually plus motion etc.)
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % units in which to provide timing for the parts of the design matrix that are specified manually
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.51; % TR
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 48; % number of slices
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 48/2; % align predictor variables so they are predicting responses midway through volume acquisition
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/mni/',subcode{s},'/'],['s6',subcode{s},'_task-namingslowish_acq-MEMB_rec-tedana_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {[root,'/SPM/first/mni/',subcode{s},'/motion_run-01.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-01_type-univariate_design-matrix.mat']}; % design matrix specific to this run
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128; % high-pass filter - slow signal drifts with a period LONGER than this will be removed (this is the default value)
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {}); % format to expect design matrix in
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {}); % format to expect motion regressors in
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/mni/',subcode{s},'/'],['s6',subcode{s},'_task-namingslowish_acq-MEMB_rec-tedana_run-02_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi_reg = {[root,'/SPM/first/mni/',subcode{s},'/motion_run-02.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-02_type-univariate_design-matrix.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).hpf = 128;
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/mni/',subcode{s},'/'],['s6',subcode{s},'_task-namingslowish_acq-MEMB_rec-tedana_run-03_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).multi_reg = {[root,'/SPM/first/mni/',subcode{s},'/motion_run-03.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-03_type-univariate_design-matrix.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).hpf = 128; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(3).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).scans = cellstr(spm_select('ExtFPList',[root,'/SPM/first/mni/',subcode{s},'/'],['s6',subcode{s},'_task-namingslowish_acq-MEMB_rec-tedana_run-04_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'],[1:551])); 
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).multi_reg = {[root,'/SPM/first/mni/',subcode{s},'/motion_run-04.mat']};
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).multi = {['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode{s},'_run-04_type-univariate_design-matrix.mat']}; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).hpf = 128; 
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(4).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.mask = {'/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'}; % MNI space brain mask
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; % do not assume that the HRF varies between subjects and voxels
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1; % use linear convolution for the HRF (https://pubmed.ncbi.nlm.nih.gov/9438436/)
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None'; % no global intensity normalisation
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0; % default to go with the above
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; % model autocorrelation

            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; % as opposed to Bayesian

            matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'Task>Fix';
            matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = [[1 -1; -1 1] zeros(2,8) [1 -1; -1 1] zeros(2,8) [1 -1; -1 1] zeros(2,8) [1 -1; -1 1] zeros(2,8)];
            matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'A';
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 zeros(1,9) 1 zeros(1,9) 1 zeros(1,9) 1 zeros(1,9)];
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'I';
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 1 zeros(1,9) 1 zeros(1,9) 1 zeros(1,9) 1 zeros(1,8)];
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'A>I';
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [1 -1 zeros(1,8) 1 -1 zeros(1,8) 1 -1 zeros(1,8) 1 -1 zeros(1,8)];
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'I>A';
            matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [-1 1 zeros(1,8) -1 1 zeros(1,8) -1 1 zeros(1,8) -1 1 zeros(1,8) ];
            matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'A+I';
            matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [1 1 zeros(1,8) 1 1 zeros(1,8) 1 1 zeros(1,8) 1 1 zeros(1,8) ];
            matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
            spm_jobman('run',matlabbatch);

        end
    end
end