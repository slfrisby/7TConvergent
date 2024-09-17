#!/bin/bash

export PATH=/imaging/local/software/centos7/ants/bin/ants/bin/:$PATH
export ANTSPATH=/imaging/local/software/centos7/ants/bin/ants/bin/

subcode=$1

dirp=/imaging/projects/cbu/wbic-p00591-DAISY/main
work=/imaging/projects/cbu/wbic-p00591-DAISY/main/work

antsApplyTransforms -d 3 -i $work/$subcode/cox_coord.nii -r $dirp/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii -o $work/$subcode/cox_coord_transformed.nii --default-value 0 --float 1 -n NearestNeighbor --transform $dirp/derivatives/fmriprep/$subcode/anat/"$subcode"_run-1_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5 --transform [$dirp/derivatives/fmriprep/$subcode/anat/"$subcode"_run-1_from-T1w_to-native_mode-image_xfm.mat,1] --transform identity --transform identity






