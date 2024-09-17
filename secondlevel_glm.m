addpath('/group/mlr-lab/AH/Projects/spm12/');
addpath('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts');

root='/imaging/projects/cbu/wbic-p00591-DAISY/main/derivatives';
cd(root);

% make output directories
mkdir([root,'/GLM/second/maineffect/t2star/'])
mkdir([root,'/GLM/second/maineffect/tedana/'])
mkdir([root,'/GLM/second/animate/t2star/'])
mkdir([root,'/GLM/second/animate/tedana/'])
mkdir([root,'/GLM/second/inanimate/t2star/'])
mkdir([root,'/GLM/second/inanimate/tedana/'])
mkdir([root,'/GLM/second/animate_vs_inanimate/t2star/'])
mkdir([root,'/GLM/second/animate_vs_inanimate/tedana/'])

% good subjects only
subcode = {'sub-001','sub-002','sub-003','sub-004','sub-007','sub-009','sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019','sub-020','sub-021','sub-022','sub-023','sub-024','sub-026','sub-028','sub-029','sub-030','sub-031','sub-032'};

% t2star

% one-sample t-test - main effect (A+I)

clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = {[root,'/GLM/second/maineffect/t2star']};
% load animate > inanimate contrast for all participants
for s = 1:length(subcode)
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{s,1} = [root,'/GLM/first/mni/t2star/',subcode{s},'/con_0006.nii,1'];
end
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {}); % format to expect covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {}); % format to expect multiple covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1; % do not use threshold masking 
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1; % use implicit mask (exclude NaNs)
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1; % options for PET and VBM (not fMRI)
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; % as opposed to Bayesian

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'A+I';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
spm_jobman('run',matlabbatch);

% one-sample t-test - animate only

clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = {[root,'/GLM/second/animate/t2star']};
% load animate > inanimate contrast for all participants
for s = 1:length(subcode)
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{s,1} = [root,'/GLM/first/mni/t2star/',subcode{s},'/con_0002.nii,1'];
end
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {}); % format to expect covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {}); % format to expect multiple covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1; % do not use threshold masking 
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1; % use implicit mask (exclude NaNs)
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1; % options for PET and VBM (not fMRI)
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; % as opposed to Bayesian

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Animate';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
spm_jobman('run',matlabbatch);

% one-sample t-test - inanimate only

clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = {[root,'/GLM/second/inanimate/t2star']};
% load animate > inanimate contrast for all participants
for s = 1:length(subcode)
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{s,1} = [root,'/GLM/first/mni/t2star/',subcode{s},'/con_0003.nii,1'];
end
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {}); % format to expect covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {}); % format to expect multiple covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1; % do not use threshold masking 
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1; % use implicit mask (exclude NaNs)
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1; % options for PET and VBM (not fMRI)
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; % as opposed to Bayesian

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Inanimate';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
spm_jobman('run',matlabbatch);

% one-sample t-test - animate vs. inanimate

clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = {[root,'/GLM/second/animate_vs_inanimate/t2star']};
% load animate > inanimate contrast for all participants
for s = 1:length(subcode)
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{s,1} = [root,'/GLM/first/mni/t2star/',subcode{s},'/con_0004.nii,1'];
end
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {}); % format to expect covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {}); % format to expect multiple covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1; % do not use threshold masking 
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1; % use implicit mask (exclude NaNs)
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1; % options for PET and VBM (not fMRI)
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; % as opposed to Bayesian

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'A>I';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'I>A';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
spm_jobman('run',matlabbatch);

% tedana

% one-sample t-test - main effect (A+I)

clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = {[root,'/GLM/second/maineffect/tedana']};
% load animate > inanimate contrast for all participants
for s = 1:length(subcode)
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{s,1} = [root,'/GLM/first/mni/tedana/',subcode{s},'/con_0006.nii,1'];
end
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {}); % format to expect covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {}); % format to expect multiple covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1; % do not use threshold masking 
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1; % use implicit mask (exclude NaNs)
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1; % options for PET and VBM (not fMRI)
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; % as opposed to Bayesian

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'A+I';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
spm_jobman('run',matlabbatch);

% one-sample t-test - animate only

clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = {[root,'/GLM/second/animate/tedana']};
% load animate > inanimate contrast for all participants
for s = 1:length(subcode)
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{s,1} = [root,'/GLM/first/mni/tedana/',subcode{s},'/con_0002.nii,1'];
end
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {}); % format to expect covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {}); % format to expect multiple covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1; % do not use threshold masking 
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1; % use implicit mask (exclude NaNs)
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1; % options for PET and VBM (not fMRI)
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; % as opposed to Bayesian

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Animate';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
spm_jobman('run',matlabbatch);

% one-sample t-test - inanimate only

clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = {[root,'/GLM/second/inanimate/tedana']};
% load animate > inanimate contrast for all participants
for s = 1:length(subcode)
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{s,1} = [root,'/GLM/first/mni/tedana/',subcode{s},'/con_0003.nii,1'];
end
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {}); % format to expect covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {}); % format to expect multiple covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1; % do not use threshold masking 
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1; % use implicit mask (exclude NaNs)
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1; % options for PET and VBM (not fMRI)
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; % as opposed to Bayesian

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Inanimate';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
spm_jobman('run',matlabbatch);

% one-sample t-test - animate vs. inanimate

clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = {[root,'/GLM/second/animate_vs_inanimate/tedana']};
% load animate > inanimate contrast for all participants
for s = 1:length(subcode)
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{s,1} = [root,'/GLM/first/mni/tedana/',subcode{s},'/con_0004.nii,1'];
end
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {}); % format to expect covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {}); % format to expect multiple covariates in (there are none)
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1; % do not use threshold masking 
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1; % use implicit mask (exclude NaNs)
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1; % options for PET and VBM (not fMRI)
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; % as opposed to Bayesian

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'A>I';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'I>A';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
spm_jobman('run',matlabbatch);


