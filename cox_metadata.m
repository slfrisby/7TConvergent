% by Saskia. Makes metadata for all 7T decoding workflows.
% 31.05.24 - also makes PERMUTATION_STRUCT for WISC MVPA. 

% At the bottom of this script is code for calculating the subdivision of
% the ATL ROI into anterior and posterior sections.

root = ['/imaging/projects/cbu/wbic-p00591-DAISY/main/derivatives'];
cd([root]);

subcode = {'sub-001','sub-002','sub-003','sub-004','sub-007','sub-009','sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019','sub-020','sub-021','sub-022','sub-023','sub-024','sub-026','sub-028','sub-029','sub-030','sub-031','sub-032'};
data = {'t2star','tedana'};

load('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/stimulimaster.mat')
load('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/dilkina_norms.mat')
load('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/cvind.mat')

% for each data type
for d = 1:length(data)

    %  make full metadata (400 rows). Initialise output
    metadata = {};

    % for every participant
    for s = 1:length(subcode)

        % load data matrix (X)
        load([root,'/cox/',subcode{s},'/',subcode{s},'_rec-',data{d},'_X.mat']);
        % load coordinates
        load(['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/',subcode{s},'/cox_coords_',data{d},'.mat'])
        % load ROI filters
        load(['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/',subcode{s},'/L_ROI_idx_',data{d},'.mat'])
        load(['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/',subcode{s},'/R_ROI_idx_',data{d},'.mat'])
        
        % subject ID
        metadata(s).subject = str2num(erase(subcode{s},'sub-'));
        % targets - animate/inanimate. N.B. ANIMATE = 0, INANIMATE = 1 (to
        % match ECoG) 
        metadata(s).targets(1).label = 'animacy';
        metadata(s).targets(1).type = 'category';
        metadata(s).targets(1).sim_source = [];
        metadata(s).targets(1).sim_metric = [];
        metadata(s).targets(1).target = repmat([zeros(50,1);ones(50,1)],4,1);
        % targets - Dilkina norms
        metadata(s).targets(2).label = 'semantic';
        metadata(s).targets(2).type = 'similarity';
        metadata(s).targets(2).sim_source = 'Dilkina_Normalized';
        metadata(s).targets(2).sim_metric = 'cosine';
        metadata(s).targets(2).target = dilkina_norms; % N.B. these are not extended to 400 stimuli
        % stimuli
        metadata(s).stimuli = repmat(stimulimaster,4,1);
        % filters
        % row and column filters currently include all rows or columns
        metadata(s).filters(1).label = 'rowfilter';
        metadata(s).filters(1).dimension = 1;
        metadata(s).filters(1).filter = logical(ones(400,1));
        metadata(s).filters(2).label = 'colfilter';
        metadata(s).filters(2).dimension = 2;
        metadata(s).filters(2).filter = logical(ones(1,size(X,2)));
        % animate and inanimate row filters
        metadata(s).filters(3).label = 'animate';
        metadata(s).filters(3).dimension = 1;
        metadata(s).filters(3).filter = repmat(logical([ones(50,1);zeros(50,1)]),4,1);
        metadata(s).filters(4).label = 'inanimate';
        metadata(s).filters(4).dimension = 1;
        metadata(s).filters(4).filter = repmat(logical([zeros(50,1);ones(50,1)]),4,1);
        % ROI filters - whole ROIs 
        metadata(s).filters(5).label = 'left_ATL_ROI';
        metadata(s).filters(5).dimension = 2;
        metadata(s).filters(5).filter = logical(L_ROI_idx);
        metadata(s).filters(6).label = 'right_ATL_ROI';
        metadata(s).filters(6).dimension = 2;
        metadata(s).filters(6).filter = logical(R_ROI_idx);
        % based on the above (and on maths at the bottom of this script),
        % create anterior and posterior halves of the ROI
        % initialise
        L_ROI_idx_ant = L_ROI_idx;
        L_ROI_idx_pos = L_ROI_idx;
        R_ROI_idx_ant = R_ROI_idx;
        R_ROI_idx_pos = R_ROI_idx;
        % for each trio of coordinates
        for i = 1:size(coords,1)
            % if the y-coordinate is greater than or equal to -20
            if coords(i,2) >= -20
                % it isn't in the posterior half of the ROI
                L_ROI_idx_pos(i) = 0;
                R_ROI_idx_pos(i) = 0;
            % if the y-coordiate is less than -20
            elseif coords(i,2) < -20
                % it isn't in the anterior half
                L_ROI_idx_ant(i) = 0;
                R_ROI_idx_ant(i) = 0;
            end
        end
        metadata(s).filters(7).label = 'left_ATL_ROI_ant';
        metadata(s).filters(7).dimension = 2;
        metadata(s).filters(7).filter = logical(L_ROI_idx_ant);
        metadata(s).filters(8).label = 'left_ATL_ROI_pos';
        metadata(s).filters(8).dimension = 2;
        metadata(s).filters(8).filter = logical(L_ROI_idx_pos);
        metadata(s).filters(9).label = 'right_ATL_ROI_ant';
        metadata(s).filters(9).dimension = 2;
        metadata(s).filters(9).filter = logical(R_ROI_idx_ant);
        metadata(s).filters(10).label = 'right_ATL_ROI_pos';
        metadata(s).filters(10).dimension = 2;
        metadata(s).filters(10).filter = logical(R_ROI_idx_pos);
        % coordinates
        metadata(s).coords.orientation = 'mni';
        metadata(s).coords.labels = [];
        metadata(s).coords.ijk = [];
        metadata(s).coords.ind = [];
        metadata(s).coords.xyz = coords;
        % cross-validation indices
        metadata(s).cvind = repmat(cvind,4,1);
        % rows and columns
        [metadata(s).nrow, metadata(s).ncol] = size(X);
        % sessions
        metadata(s).sessions = [ones(100,1);2*ones(100,1);3*ones(100,1);4*ones(100,1);];
    end
    save([root,'/cox/',data{d},'_full_metadata.mat'],'metadata','-v7.3')

    % make metadata and PERMUTATION_STRUCT for data averaged across runs
    metadata = {};
    PERMUTATION_STRUCT = {};

    % for every participant
    for s = 1:length(subcode)

        % load data matrix (X)
        load([root,'/cox/',subcode{s},'/',subcode{s},'_rec-',data{d},'_X.mat']);
        % load coordinates
        load(['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/',subcode{s},'/cox_coords_',data{d},'.mat'])
        % load ROI filters
        load(['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/',subcode{s},'/L_ROI_idx_',data{d},'.mat'])
        load(['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/',subcode{s},'/R_ROI_idx_',data{d},'.mat'])
        
        % subject ID
        metadata(s).subject = str2num(erase(subcode{s},'sub-'));
        % targets - animate/inanimate. N.B. ANIMATE = 0, INANIMATE = 1 (to
        % match ECoG) 
        metadata(s).targets(1).label = 'animacy';
        metadata(s).targets(1).type = 'category';
        metadata(s).targets(1).sim_source = [];
        metadata(s).targets(1).sim_metric = [];
        metadata(s).targets(1).target = [zeros(50,1);ones(50,1)];
        % targets - Dilkina norms
        metadata(s).targets(2).label = 'semantic';
        metadata(s).targets(2).type = 'similarity';
        metadata(s).targets(2).sim_source = 'Dilkina_Normalized';
        metadata(s).targets(2).sim_metric = 'cosine';
        metadata(s).targets(2).target = dilkina_norms; 
        % stimuli
        metadata(s).stimuli = stimulimaster;
        % filters
        % row and column filters currently include all rows or columns
        metadata(s).filters(1).label = 'rowfilter';
        metadata(s).filters(1).dimension = 1;
        metadata(s).filters(1).filter = logical(ones(100,1));
        metadata(s).filters(2).label = 'colfilter';
        metadata(s).filters(2).dimension = 2;
        metadata(s).filters(2).filter = logical(ones(1,size(X,2)));
        % animate and inanimate row filters
        metadata(s).filters(3).label = 'animate';
        metadata(s).filters(3).dimension = 1;
        metadata(s).filters(3).filter = logical([ones(50,1);zeros(50,1)]);
        metadata(s).filters(4).label = 'inanimate';
        metadata(s).filters(4).dimension = 1;
        metadata(s).filters(4).filter = logical([zeros(50,1);ones(50,1)]);
        % ROI filters
        metadata(s).filters(5).label = 'left_ATL_ROI';
        metadata(s).filters(5).dimension = 2;
        metadata(s).filters(5).filter = logical(L_ROI_idx);
        metadata(s).filters(6).label = 'right_ATL_ROI';
        metadata(s).filters(6).dimension = 2;
        metadata(s).filters(6).filter = logical(R_ROI_idx);
        % based on the above (and on maths at the bottom of this script),
        % create anterior and posterior halves of the ROI
        % initialise
        L_ROI_idx_ant = L_ROI_idx;
        L_ROI_idx_pos = L_ROI_idx;
        R_ROI_idx_ant = R_ROI_idx;
        R_ROI_idx_pos = R_ROI_idx;
        % for each trio of coordinates
        for i = 1:size(coords,1)
            % if the y-coordinate is greater than or equal to -20
            if coords(i,2) >= -20
                % it isn't in the posterior half of the ROI
                L_ROI_idx_pos(i) = 0;
                R_ROI_idx_pos(i) = 0;
            % if the y-coordiate is less than -20
            elseif coords(i,2) < -20
                % it isn't in the anterior half
                L_ROI_idx_ant(i) = 0;
                R_ROI_idx_ant(i) = 0;
            end
        end
        metadata(s).filters(7).label = 'left_ATL_ROI_ant';
        metadata(s).filters(7).dimension = 2;
        metadata(s).filters(7).filter = logical(L_ROI_idx_ant);
        metadata(s).filters(8).label = 'left_ATL_ROI_pos';
        metadata(s).filters(8).dimension = 2;
        metadata(s).filters(8).filter = logical(L_ROI_idx_pos);
        metadata(s).filters(9).label = 'right_ATL_ROI_ant';
        metadata(s).filters(9).dimension = 2;
        metadata(s).filters(9).filter = logical(R_ROI_idx_ant);
        metadata(s).filters(10).label = 'right_ATL_ROI_pos';
        metadata(s).filters(10).dimension = 2;
        metadata(s).filters(10).filter = logical(R_ROI_idx_pos);
        % coordinates
        metadata(s).coords.orientation = 'mni';
        metadata(s).coords.labels = [];
        metadata(s).coords.ijk = [];
        metadata(s).coords.ind = [];
        metadata(s).coords.xyz = coords;
        % cross-validation indices
        metadata(s).cvind = cvind;
        % rows and columns
        metadata(s).nrow = 100;
        metadata(s).ncol = size(X,2);
        % sessions
        metadata(s).sessions = ones(100,1);

        % subject ID
        PERMUTATION_STRUCT(s).subject = str2num(erase(subcode{s},'sub-'));
        % stimuli
        PERMUTATION_STRUCT(s).stimuli.stimulus = stimulimaster;
        % permutation order 
        permutation_index = zeros(size(stimulimaster,1),10000);
        for i = 1:10000
            permutation_index(:,i) = randperm(size(stimulimaster,1))';
        end
        PERMUTATION_STRUCT(s).permutation_index = permutation_index;

    end
    save([root,'/cox/',data{d},'_averaged_metadata.mat'],'metadata','-v7.3')
    save([root,'/cox/',data{d},'_averaged_PERMUTATION_STRUCT.mat'],'PERMUTATION_STRUCT','-v7.3')

 
end

%% ATL subdivision

% % add path
% addpath('/group/mlr-lab/AH/Projects/spm12/')

% % load ROI
% vol = spm_vol('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/ROI/L_vATL_ROI.nii');
% ROI = spm_read_vols(vol);

% % load grey matter mask in MNI space. N.B. all coordinates of voxels are
% % given in MNI space, even though the voxel values are in native space.
% vol = spm_vol('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_class-GM_probtissue.nii');
% mask = spm_read_vols(vol);

% % binarise above threshold
% mask(mask > 0.2) = 1;
% mask(mask < 0.2) = 0;

% % find intersection between ROI and mask
% tmp = mask + ROI;
% tmp(tmp < 2) = 0;

% % find maximum and minimum y-values (i.e. furthest forward and furthest
% % back points of ROI)
% [x,y,z] = ind2sub(size(tmp),find(tmp~=0));
% miny = min(y);
% maxy = max(y);

% % find halfway point between them
% line = miny + 0.5*(maxy-miny);

% % convert to MNI space. Select arbitrary point where y = line
% idx = find(y==113,1);
% c = cor2mni([x(idx),y(idx),z(idx)],vol.mat);

% c(2)

% % The line bisecting the ATL ROI in the y-direction is at y = -20 in MNI
% % space.
