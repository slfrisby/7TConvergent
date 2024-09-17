function cox_format(ids)

% by Saskia. Reformats data for use in decoding, including in CHTC
% workflows. 

    addpath('/group/mlr-lab/AH/Projects/spm12/');
    addpath('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts');

    root = ['/imaging/projects/cbu/wbic-p00591-DAISY/main/derivatives'];
    cd([root]);

    subcode{1} = ['sub-',ids];
    data = {'t2star','tedana'};
    runs = {'01','02','03','04'};
    load('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/stimulimaster.mat');
    
    % for every participant
    for s = 1:size(subcode,2)

        % make output directory
        mkdir([root,'/cox/',subcode{s},'/']);

        % load stimulus orders and accuracies
        excelfile = readtable('/imaging/projects/cbu/wbic-p00591-DAISY/main/behavioural/accuracies.xlsx','Sheet',ids);

        for d = 1:length(data)
            % initialise output - (100 stimuli x 4 runs) x (84 x 84 x 48
            % voxels)
            X = zeros(400,338688);

            % make list of coordinates in matrix space. After the beta image is
            % flattened, this is the order that the voxels will be in.
            x = repmat([1:84]',84*48,1);
            y = repmat(repelem(1:84,84)',48,1);
            z = repelem([1:48]',84*84);
            coords = [x,y,z];
            clear x y z

    
            % for every run
            for r = 1:length(runs)
    
                % get order of stimuli within that run
                stimuli = table2array(excelfile(1:100,3*r-2));
                stimuli = erase(stimuli,'StimFiles/');
                stimuli = erase(stimuli,'.bmp');
    
                % also get accuracies
                accuracies = table2array(excelfile(1:100,3*r-1));
    
                % for each stimulus
                for stim = 1:100
    
                   % find the index of that stimulus in the stimulimaster file
                   idx = find(strcmp(stimulimaster,stimuli{stim}));
        
                   % if the stimulus was named incorrectly, do not load beta
                   % image. Instead, fill in the row with NaNs. These will be
                   % median-interpolated later. 
                   if ~accuracies(stim)
                       X(idx+100*(r-1),:) = NaN;
    
                   else
                        % load beta image. Note that, after every 6 stimulus betas,
                        % there are 6 motion betas, so we need to skip over those.
                        vol = spm_vol(sprintf([root,'/GLM/first/native/',data{d},'/',subcode{s},'/beta_%04d.nii'],stim+106*(r-1)));
                        img = spm_read_vols(vol);
        
                        % flatten - to undo you should use rrimg = reshape(rimg,84,84,48);
                        rimg = reshape(img,1,[]);
    
                        % put the betas into the output matrix in the order
                        % specified in stimulimaster
                        X(idx+100*(r-1),:) = rimg;
                   end
                end
            end
    
            % get just voxels within grey matter. Those outside the grey
            % matter will have NaN for every trial 
            nanidx = [];
            for i = 1:size(X,2)
                % if all the numbers in the column are NaN
                if sum(isnan(X(:,i))) == 400
                    % flag the column
                    nanidx = [nanidx i];
                end
            end
    
            % index
            X(:,nanidx) = [];
            % also index coordinates
            coords(nanidx,:) = [];
    
            % median-interpolate trials which participants got wrong. Use
            % median across all trials for each voxel to avoid incorporating
            % category information into the data (which we would do if we took
            % the median within condition)
            tmp = find(isnan(X));
            [stim, ~] = ind2sub(size(X),tmp);
            stim = unique(stim);
            interp = nanmedian(X);
            for i = 1:length(stim)
                X(stim(i),:) = interp;
            end
        
            % save data
            save([root,'/cox/',subcode{s},'/',subcode{s},'_rec-',data{d},'_X.mat'],'X','-v7.3')
            
    
            % create index of which voxels are in ROI. This will become part of
            % the metadata. 
            % 
            % % load left ROI
            % vol = spm_vol([root,'/fmriprep/',subcode{s},'/anat/',subcode{s},'_run-1_space-native_label-LROI_mask.nii.gz']);
            % L_ROI = spm_read_vols(vol);
            % L_ROI_idx = zeros(1,size(X,2));
            % % load right ROI
            % vol = spm_vol([root,'/fmriprep/',subcode{s},'/anat/',subcode{s},'_run-1_space-native_label-RROI_mask.nii.gz']);
            % R_ROI = spm_read_vols(vol);
            % R_ROI_idx = zeros(1,size(X,2));
            % for c = 1:size(coords,1)
            %     % if coordinates lie within either ROI, index 1 for that ROI
            %     if L_ROI(coords(c,1),coords(c,2),coords(c,3))
            %         L_ROI_idx(c) = 1;
            %     end
            %     if R_ROI(coords(c,1),coords(c,2),coords(c,3))
            %         R_ROI_idx(c) = 1;
            %     end
            % end
            % % save
            % save(['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/',subcode{s},'/L_ROI_idx_',data{d},'.mat'],'L_ROI_idx','-v7.3');
            % save(['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/',subcode{s},'/R_ROI_idx_',data{d},'.mat'],'R_ROI_idx','-v7.3');
            % 
            % SOSLASSO takes data in matrix space, but coordinates that give
            % the corresponding location of each voxel in MNI space.
            
            % % load example functional volume to obtain header information
            % vol = spm_vol([root,'/halaiprep/',subcode{s},'/func/',subcode{s},'_task-namingslowish_acq-MEMB_rec-',data{d},'_run-01_space-native_desc-preproc_bold.nii.gz']);
            % vol = vol(1);
            % % for every set of coordinates in the grey matter
            % for c = 1:size(coords,1)
            %     % construct an empty volume with just that voxel filled in
            %     img = zeros(84,84,48);
            %     img(coords(c,1),coords(c,2),coords(c,3)) = 1;
            %     % write volume 
            %     tmp = vol;
            %     tmp.fname = ['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/',subcode{s},'/cox_coord.nii'];
            %     spm_write_vol(tmp,img);
            %     % call bash script that uses AntsApplyTransforms to transform
            %     % coordinate into MNI space
            %     system(['bash /imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/cox_find_coords.sh ',subcode{s}]);
            %     % read transformed file back in 
            %     tmp = spm_vol(['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/',subcode{s},'/cox_coord_transformed.nii']);
            %     img = spm_read_vols(tmp);
            %     % get T matrix
            %     T = tmp.mat;
            %     % find 1s in volume (i.e. where the voxel has gone to)
            %     [x,y,z] = ind2sub(size(img),find(img==1));
            %     % to find the centre of the sphere, take the (rounded) mean x,
            %     % y, and z coordinates. Run these through cor2mni and overwrite
            %     % the matrix space coordinates
            %     coords(c,:) = cor2mni([round(mean(x)), round(mean(y)), round(mean(z))],T);
            % end
            % save
            % save(['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/',subcode{s},'/cox_coords_',data{d},'.mat'],'coords','-v7.3');
        end
    end
end

