addpath('/group/mlr-lab/AH/Projects/spm12/')
addpath('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts')
mkdir('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/ROI/')

% load template in MNI space
vol = spm_vol('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_T1w_brain.nii');
template = spm_read_vols(vol);
% save the transformation matrix
T = vol.mat;

% some of the electrode grids wrap up around the frontal lobe. We only want
% to include temporal lobe electrodes in the ROI. So load atlas
tmp = spm_vol('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/Comb_atlas.nii');
atlas = spm_read_vols(tmp);
% also save its transformation matrix
atlasT = tmp.mat;

% construct empty image in which to construct ROI
roi = nan(size(template));

% load electrode locations 
coords = readtable('MNI_basal_electrodes_sub-01-22_w_label.csv');
coords = table2array(coords(:,3:5));
% remove nans
coords = rmmissing(coords,1);
% check whether electrodes are in temporal lobe. Transform coordinates into
% atlas matrix space 
atlascoords = mni2cor(coords,atlasT);
% there isn't much white space around the atlas, so very anterior/ventral
% electrodes can end up with coordinates of zero. If this is the case,
% simply fill in 1 (I have verified via inspection that this gives a
% sensible result)
atlascoords(find(atlascoords==0)) = 1;
% the same is true for negative values
atlascoords(find(atlascoords<0)) = 1;
% read atlas
atlasvals = zeros(size(coords,1),1);
for i = 1:size(coords,1)
    atlasvals(i) = atlas(atlascoords(i,1),atlascoords(i,2),atlascoords(i,3));
end
% eliminate electrodes that lie outside the temporal lobe. Notes:
% - some electrodes have values of 0 (outside brain) on the atlas - these
% have all been visualised and denote an electrode on the temporal lobe
% - one electrode appears to be in the brainstem, but is probably in the
% parahippocampal gyrus and so is included
% - electrodes in temporo-occipital regions (as labelled by the atlas) are
% included 
% - electrodes in the precentral gyrus, postcentral gyrus, supramarginal
% gyrus, and occipital fusiform are excluded
coords(find(atlasvals == 7 | atlasvals == 17 | atlasvals == 19 | atlasvals == 40),:) = [];
% convert to matrix coordinates
coords = mni2cor(coords,T);

% loop electrodes. set voxels where electrodes are located to a high number
% to ensure that the spatial extent of the smoothing is preserved
for elec = 1:size(coords,1)
    roi(coords(elec,1),coords(elec,2),coords(elec,3)) = 1;
end

% save unsmoothed image. Using niftiwrite preserves nans in the image,
% which is important - using zeros will mean that the sparse 1s used to
% indicate electrode location are destroyed during smoothing. The header
% information will no longer be correct, so later on we force the image to
% be saved with the same header info as the template.
niftiwrite(roi,'/imaging/projects/cbu/wbic-p00591-DAISY/main/work/electrode_locations.nii');

% smooth - 10mm FWHM creates a smooth shape (this figure decided by trial
% and error)
clear matlabbatch
matlabbatch{1}.spm.spatial.smooth.data = cellstr(['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/electrode_locations.nii']);
matlabbatch{1}.spm.spatial.smooth.fwhm = [10 10 10];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 'smoothed_';
spm_jobman('run',matlabbatch);

% load smoothed file
tmp = spm_vol('/imaging/projects/cbu/wbic-p00591-DAISY/main/work/smoothed_electrode_locations.nii');
roi = spm_read_vols(tmp);

% binarise above liberal threshold
roi(roi > 0.0001) = 1;

% get left ROI only
leftroi = roi;
% set all right hemisphere values to zero
leftroi(96:193,:,:)=0;
% Generous smoothing means that we have wiggle room
% around the temporal lobe, which will be useful in backprojection.
% However, we do not want to include regions outside the temporal lobe.
% create atlas mask
atlasmask = atlas;
% set left temporal lobe regions to zero
regionidx = [8,9,10,11,12,13,14,15,16,34,35,37,38,39];
atlasmask(ismember(atlasmask,regionidx)) = 0;
% set other regions to 1
atlasmask(atlasmask>0) = 1;
% the Harvard-Oxford atlas does not include a cerebellum, but is overall
% useful because it is big enough to label electrodes on the surface of the
% brain accurately. The AAL atlas is slightly smaller than the MNI brain,
% meaning that it is difficult to use to label electrodes. However, it does
% have a cerebellum. So load the AAL atlas, get the cerebellar regions, and
% add them to the mask
tmp = spm_vol('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/waal_MNI_V4.nii');
aal = spm_read_vols(tmp);
cerebellumidx = [91:108];
atlasmask(ismember(aal,cerebellumidx)) = 1;

for x = 1:size(leftroi,1)
    for y = 1:size(leftroi,2)
        for z = 1:size(leftroi,3)
            % if the voxel is in the ROI
            if leftroi(x,y,z) == 1
                % translate the coordinates into atlas space
                ac = mni2cor(cor2mni([x,y,z],T),atlasT);
                % since the template has more whitespace around it than the
                % atlas, we sometimes get ac values that are larger or smaller than
                % the atlas matrix. Force these to be within the atlaas
                % matrix
                if ac(1) > size(atlas,1)
                    ac(1) = size(atlas,1);
                elseif ac(1) < 1
                    ac(1) = 1;
                end
                if ac(2) > size(atlas,2)
                    ac(2) = size(atlas,2);
                elseif ac(2) < 1
                    ac(2) = 1;
                end
                if ac(3) > size(atlas,3)
                    ac(3) = size(atlas,3);
                elseif ac(3) < 1
                    ac(3) = 1;
                end
                % if the atlas coordinates lie within the atlas mask
                if atlasmask(ac(1),ac(2),ac(3)) == 1
                    % set the ROI at that voxel to zero
                    leftroi(x,y,z) = 0;
                end
            end
        end
    end
end

% save temporarily
tmp = vol;
tmp.fname = '/imaging/projects/cbu/wbic-p00591-DAISY/main/work/tmp.nii';
spm_write_vol(tmp,leftroi);
gzip(tmp.fname);

% smooth again to remove jagged edges and small fragments associated
% with masking (this figure decided by trial and error)
clear matlabbatch
matlabbatch{1}.spm.spatial.smooth.data = cellstr(['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/tmp.nii']);
matlabbatch{1}.spm.spatial.smooth.fwhm = [4 4 4];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 'smoothed_';
spm_jobman('run',matlabbatch);

% load file
tmp = spm_vol('/imaging/projects/cbu/wbic-p00591-DAISY/main/work/smoothed_tmp.nii');
leftroi = spm_read_vols(tmp);

% write left ROI
tmp = vol;
tmp.fname = '/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/ROI/L_vATL_ROI.nii';
spm_write_vol(tmp,leftroi);
gzip(tmp.fname);

% because there are fewer right hemisphere patients, the right ROI is
% smaller. We want the same coverage in both hemispheres, so flip the left
% ROI into the right hemisphere
rightroi = leftroi(193:-1:1,:,:);
% save - force header information to match MNI template
tmp = vol;
tmp.fname = '/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/ROI/R_vATL_ROI.nii';
spm_write_vol(tmp,rightroi);
gzip(tmp.fname);

% for ease of checking results, create an image of the electrode
% locations, smoothed at 2 mm FWHM for visibility (and because the
% recording diameter of the electrodes is 2.3 mm). The smoothed file from
% earlier is no longer needed and can be overwritten.
clear matlabbatch
matlabbatch{1}.spm.spatial.smooth.data = cellstr(['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/electrode_locations.nii']);
matlabbatch{1}.spm.spatial.smooth.fwhm = [2 2 2];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 'smoothed_';
spm_jobman('run',matlabbatch);

% load electrode image
tmp = spm_vol('/imaging/projects/cbu/wbic-p00591-DAISY/main/work/smoothed_electrode_locations.nii');
elecs = spm_read_vols(tmp);
% binarise above threshold (so that SPM can save it correctly - threshold selected for clarity of resulting figure)
elecs(elecs > 0.001) = 1;
% write with corrected header
tmp = vol;
tmp.fname = '/imaging/projects/cbu/wbic-p00591-DAISY/main/work/smoothed_electrode_locations.nii';
spm_write_vol(tmp,elecs);



