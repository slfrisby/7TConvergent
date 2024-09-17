#!/bin/bash

# load workbench
module load workbench/1.5.0

# set directories
dirp=/imaging/projects/cbu/wbic-p00591-DAISY/main
work=/imaging/projects/cbu/wbic-p00591-DAISY/main/work
# $work/coefficients should already be present because output from CHTC should be saved there

# unzip template mesh from fsaverage
if [ ! -f $work/coefficients/pial_right.gii ]; then
cp $dirp/scripts/fsaverage/pial_left.gii.gz $work/coefficients/pial_left.gii.gz
cp $dirp/scripts/fsaverage/pial_right.gii.gz $work/coefficients/pial_right.gii.gz
gunzip $work/coefficients/pial_left.gii.gz
gunzip $work/coefficients/pial_right.gii.gz
fi

mkdir $work/ROIsurf/

# project to the pial surface

# Left ROI - whole, anterior, posterior
wb_command -volume-to-surface-mapping $dirp/scripts/ROI/L_vATL_ROI.nii $work/coefficients/pial_left.gii $work/ROIsurf/L_vATL_ROI.func.gii -trilinear
wb_command -volume-to-surface-mapping $dirp/scripts/ROI/L_vATL_ant.nii $work/coefficients/pial_left.gii $work/ROIsurf/L_vATL_ant_ROI.func.gii -trilinear
wb_command -volume-to-surface-mapping $dirp/scripts/ROI/L_vATL_pos.nii $work/coefficients/pial_left.gii $work/ROIsurf/L_vATL_pos_ROI.func.gii -trilinear

# Right ROI - whole, anterior, posterior
wb_command -volume-to-surface-mapping $dirp/scripts/ROI/R_vATL_ROI.nii $work/coefficients/pial_right.gii $work/ROIsurf/R_vATL_ROI.func.gii -trilinear
wb_command -volume-to-surface-mapping $dirp/scripts/ROI/R_vATL_ant.nii $work/coefficients/pial_right.gii $work/ROIsurf/R_vATL_ant_ROI.func.gii -trilinear
wb_command -volume-to-surface-mapping $dirp/scripts/ROI/R_vATL_pos.nii $work/coefficients/pial_right.gii $work/ROIsurf/R_vATL_pos_ROI.func.gii -trilinear

# Univariate contrasts
wb_command -volume-to-surface-mapping $dirp/derivatives/GLM/second/animate_vs_inanimate/tedana/A_gt_I_0001_FWE05.nii $work/coefficients/pial_left.gii $work/ROIsurf/L_univariate_animate.func.gii -trilinear
wb_command -volume-to-surface-mapping $dirp/derivatives/GLM/second/animate_vs_inanimate/tedana/A_gt_I_0001_FWE05.nii $work/coefficients/pial_right.gii $work/ROIsurf/R_univariate_animate.func.gii -trilinear
wb_command -volume-to-surface-mapping $dirp/derivatives/GLM/second/animate_vs_inanimate/tedana/I_gt_A_0001_FWE05.nii $work/coefficients/pial_left.gii $work/ROIsurf/L_univariate_inanimate.func.gii -trilinear
wb_command -volume-to-surface-mapping $dirp/derivatives/GLM/second/animate_vs_inanimate/tedana/I_gt_A_0001_FWE05.nii $work/coefficients/pial_right.gii $work/ROIsurf/R_univariate_inanimate.func.gii -trilinear

# tSNR map
# gunzip $dirp/derivatives/tSNR/group_tSNRmap.nii.gz
wb_command -volume-to-surface-mapping $dirp/derivatives/tSNR/group_tSNRmap.nii $work/coefficients/pial_left.gii $work/ROIsurf/L_tSNR.func.gii -trilinear
wb_command -volume-to-surface-mapping $dirp/derivatives/tSNR/group_tSNRmap.nii $work/coefficients/pial_right.gii $work/ROIsurf/R_tSNR.func.gii -trilinear

 

