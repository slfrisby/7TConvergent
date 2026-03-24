% Compare tSNR and functional contrast-to-noise between anterior and
% posterior ROIs, using data from previous 7T methods study. 

addpath('/group/mlr-lab/AH/Projects/spm12/');
addpath('/group/mlr-lab/AH/Projects/toolboxes/');
addpath('/group/mlr-lab/AH/Projects/toolboxes/Violinplot/');
addpath('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/riksneurotools-master/Util');
% the above line adds the function roi_extract.m to path. Before using this
% % function, check that lines 195-196 read:
%                 ROI(nr,s).mean   = nanmean(d,2);
%                 ROI(nr,s).median = nanmedian(d,2);
addpath('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts');

root = ['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/derivatives'];
cd(root);

% exclusions: 013 and 017 for excessive head motion, and 004 because the
% pTx run was not acquired.
subs=[{'001'},{'002'},{'003'},{'005'},{'006'},{'007'},{'008'},{'009'},{'010'},{'011'},{'012'},{'014'},{'015'},{'016'},{'017'},{'018'},{'019'}];

% setup ROIs. Use all ROIs that overlapped with any main effect mentioned
% in the main text
R = struct();
R.ROIfiles{1}=['/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/ROI/L_vATL_ant.nii']; % left ATL anterior
R.ROIfiles{2}=['/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/ROI/L_vATL_pos.nii']; % left vATL posterior
R.ROIfiles{3}=['/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/ROI/R_vATL_ant.nii'];% right vATL anterior
R.ROIfiles{4}=['/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/ROI/R_vATL_pos.nii'];% right vATL posterior

% extract ROIs from con images (expressing model fit)

imgs = cell(1,length(subs));
con = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB/con_0003.nii']; 
end

% extract data from each ROI
con.Datafiles = imgs;
con.output_raw = 1;
con.ROIfiles=R.ROIfiles;
con.ROI = roi_extract(con);

% extract ROIs from spmT images (expressing model precision)

imgs = cell(1,length(subs));
spmT = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB/spmT_0003.nii']; 
end

% extract data from each ROI
spmT.Datafiles = imgs;
spmT.output_raw = 1;
spmT.ROIfiles=R.ROIfiles;
spmT.ROI = roi_extract(spmT);

% extract ROIs from tSNR images 

imgs = cell(1,length(subs));
tSNR = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/tSNR/sub-',subs{s},'/MEMB_tSNR.nii.gz']; 
end

% extract data from each ROI
tSNR.Datafiles = imgs;
tSNR.output_raw = 1;
tSNR.ROIfiles=R.ROIfiles;
tSNR.ROI = roi_extract(tSNR);

%remove NaN voxels from data
% for each participant
for s=1:length(subs)
    % for each ROI
    for n=1:length(con.ROIfiles)
        % find NaNs in the contrast
        x=isnan(con.ROI(n,s).rawdata(1,:));
        % get indices of these NaNs
        ind=find(x);
        % remove the NaNs from the data and the corresponding values
        % from the spmT data
        con.ROI(n,s).rawdata(:,ind)=[];
        spmT.ROI(n,s).rawdata(:,ind)=[];
        tSNR.ROI(n,s).rawdata(:,ind)=[];
        % remove the coordinates of the NaNs from the matrix of
        % coordinates, so that everything lines up
        con.ROI(n,s).XYZ(:,ind)=[];
        spmT.ROI(n,s).XYZ(:,ind)=[];
        tSNR.ROI(n,s).XYZ(:,ind)=[];
        % extract the mean and median for each ROI (for each contrast)
        % !! N.B. this is done slightly differently to how it was in the
        % main methods study, in order to account for the fact that these
        % ROIs overlap non-brain regions (which skews roi_extract.m)
        con_collate_median{n}(s,:)=median(con.ROI(n,s).rawdata);
        con_collate_mean{n}(s,:)=mean(con.ROI(n,s).rawdata);
        spmT_collate_median{n}(s,:)=median(spmT.ROI(n,s).rawdata);
        spmT_collate_mean{n}(s,:)=mean(spmT.ROI(n,s).rawdata);
        tSNR_collate_median{n}(s,:)=median(tSNR.ROI(n,s).rawdata);
        tSNR_collate_mean{n}(s,:)=mean(tSNR.ROI(n,s).rawdata);
    end
end

% t-tests - fCNR/tSNR greater in posterior than anterior half? 

% contrast
% left hemisphere
% conduct paired t-test
[tmp1,tmp2,tmp3,tmp4]=ttest(con_collate_median{1},con_collate_median{2},'tail','left');
% get p-values
con_p_val(1)=tmp2;
% also get t-values
con_t_val(1)=tmp4.tstat;
% right hemisphere
% conduct paired t-test
[tmp1,tmp2,tmp3,tmp4]=ttest(con_collate_median{3},con_collate_median{4},'tail','left');
% get p-values
con_p_val(2)=tmp2;
% also get t-values
con_t_val(2)=tmp4.tstat;

% spmT
% left hemisphere
% conduct paired t-test
[tmp1,tmp2,tmp3,tmp4]=ttest(spmT_collate_median{1},spmT_collate_median{2},'tail','left');
% get p-values
spmT_p_val(1)=tmp2;
% also get t-values
spmT_t_val(1)=tmp4.tstat;
% right hemisphere
% conduct paired t-test
[tmp1,tmp2,tmp3,tmp4]=ttest(spmT_collate_median{3},spmT_collate_median{4},'tail','left');
% get p-values
spmT_p_val(2)=tmp2;
% also get t-values
spmT_t_val(2)=tmp4.tstat;

% tSNR
% left hemisphere
% conduct paired t-test
[tmp1,tmp2,tmp3,tmp4]=ttest(tSNR_collate_median{1},tSNR_collate_median{2},'tail','left');
% get p-values
tSNR_p_val(1)=tmp2;
% also get t-values
tSNR_t_val(1)=tmp4.tstat;
% right hemisphere
% conduct paired t-test
[tmp1,tmp2,tmp3,tmp4]=ttest(tSNR_collate_median{3},spmT_collate_median{4},'tail','left');
% get p-values
tSNR_p_val(2)=tmp2;
% also get t-values
tSNR_t_val(2)=tmp4.tstat;

% extract ROIs from beta images for exploratory MVPA

imgs = cell(1,length(subs));
con_mvpa = struct();

% for each participant
for s = 1:length(subs)
    % setup beta image structure
    imgs{s}{1} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0001.nii'];
    imgs{s}{2} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0002.nii'];
    imgs{s}{3} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0003.nii'];
    imgs{s}{4} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0004.nii'];
    imgs{s}{5} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0005.nii'];
    imgs{s}{6} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0006.nii'];
    imgs{s}{7} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0007.nii'];
    imgs{s}{8} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0008.nii'];
    imgs{s}{9} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0009.nii'];
    imgs{s}{10} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0010.nii'];
    imgs{s}{11} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0011.nii'];
    imgs{s}{12} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0012.nii'];
    imgs{s}{13} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0013.nii'];
    imgs{s}{14} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0014.nii'];
    imgs{s}{15} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0015.nii'];
    imgs{s}{16} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0016.nii'];
    imgs{s}{17} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0017.nii'];
    imgs{s}{18} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0018.nii'];
    imgs{s}{19} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0019.nii'];
    imgs{s}{20} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0020.nii'];
    imgs{s}{21} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0021.nii'];
    imgs{s}{22} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0022.nii'];
    imgs{s}{23} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0023.nii'];
    imgs{s}{24} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_MEMB/beta_0024.nii'];
end
    
% extract data from each ROI and store in one big struct
con_mvpa.Datafiles = imgs;
con_mvpa.output_raw = 1;
con_mvpa.ROIfiles = R.ROIfiles;
con_mvpa.ROI = roi_extract(con_mvpa);


% remove NaN voxels from data
% for each participant
for s=1:length(subs)
    % for each ROI
    for n=1:length(R.ROIfiles)
        % find NaNs in the beta images
        x=isnan(con_mvpa.ROI(n,s).rawdata(1,:));
        % get indices of these NaNs
        ind=find(x);
        % remove the NaNs from the data
        con_mvpa.ROI(n,s).rawdata(:,ind)=[];
        % remove the coordinates of the NaNs from the matrix of
        % coordinates, so that everything lines up
        con_mvpa.ROI(n,s).XYZ(:,ind)=[];
    end
end


% for each participant
for s=1:length(subs)
    % for each ROI
    for n=1:length(R.ROIfiles)
        % extract data (this matrix is beta images x nonzero voxels)
        x=con_mvpa.ROI(n,s).rawdata;
        % Each block (12 semantic and 12 control) is one row of x.
        % Calculate cosine distance between each pair of blocks (This
        % is stored as the upper triangle of the similarity matrix). 
        con_mvpa.ROI(n,s).dissimilarity=triu(squareform(pdist(x,'cosine')));
        % convert zeros in the matrix to NaNs
        con_mvpa.ROI(n,s).dissimilarity(con_mvpa.ROI(n,s).dissimilarity==0)=nan;
        % calculate the mean dissimilarity between pairs of blocks in
        % the same condition (i.e. mean dissimilarity between pairs of
        % semantic blocks or pairs of control blocks)
        con_mvpa.ROI(n,s).mvpa_within_mean=nanmean([reshape(con_mvpa.ROI(n,s).dissimilarity([1:12],[1:12]),[],1);reshape(con_mvpa.ROI(n,s).dissimilarity([13:24],[13:24]),[],1)]);
        % calculate the mean dissimilarity between pairs of blocks in
        % different conditions (i.e. mean dissimilarity between one
        % semantic block and one control block)
        con_mvpa.ROI(n,s).mvpa_between_mean=nanmean(reshape(con_mvpa.ROI(n,s).dissimilarity([1:12],[13:24]),[],1));
        % calculate the difference in means
        con_mvpa.ROI(n,s).mvpa_comparison_mean=[con_mvpa.ROI(n,s).mvpa_between_mean-con_mvpa.ROI(n,s).mvpa_within_mean];
        % collate results
        con_mvpa_collate_mean{n}(s)=con_mvpa.ROI(n,s).mvpa_comparison_mean;
    end
end

% t-test - fCNR/tSNR greater in posterior than anterior half? 

% left hemisphere
% conduct paired t-test
[tmp1,tmp2,tmp3,tmp4]=ttest(con_mvpa_collate_mean{1},con_mvpa_collate_mean{2},'tail','left');
% get p-values
con_mvpa_p_val(1)=tmp2;
% also get t-values
con_mvpa_t_val(1)=tmp4.tstat;
% right hemisphere
% conduct paired t-test
[tmp1,tmp2,tmp3,tmp4]=ttest(con_mvpa_collate_mean{3},con_mvpa_collate_mean{4},'tail','left');
% get p-values
con_mvpa_p_val(2)=tmp2;
% also get t-values
con_mvpa_t_val(2)=tmp4.tstat;

