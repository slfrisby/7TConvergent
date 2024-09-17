#!/bin/bash
echo "
++++++++++++++++++++++++" 
echo +* "Set up script run environment" 
#adds appropriate tools and options
export PATH=$PATH:/group/mlr-lab/AH/Projects/toolboxes/afni/v18.3.03
export DYLD_FALLBACK_LIBRARY_PATH=/group/mlr-lab/AH/Projects/toolboxes/afni/v18.3.03
FSLDIR=/imaging/local/software/fsl/latest/x86_64/fsl
# this line simply configures FSL
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH
dirp=/imaging/projects/cbu/wbic-p00591-DAISY/main
work=/imaging/projects/cbu/wbic-p00591-DAISY/main/work
workcond="$work"/sub-"$ids"/
cd $dirp

# make output directory
mkdir -p $dirp/derivatives/tSNR/sub-"$ids"/

# for each run

for r in 01 02 03 04; do

# calculate mean and standard deviation images (for inspection). Providing no input to 3dTstat calculates the mean
3dTstat -overwrite -prefix $dirp/derivatives/tSNR/sub-"$ids"/run-"$r"_mean.nii.gz $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_task-namingslowish_acq-MEMB_rec-t2star_run-"$r"_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz 
3dTstat -stdevNOD -overwrite -prefix $dirp/derivatives/tSNR/sub-"$ids"/run-"$r"_stdev.nii.gz $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_task-namingslowish_acq-MEMB_rec-t2star_run-"$r"_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz 
# calculate tSNR
3dTstat -tsnr -overwrite -prefix $dirp/derivatives/tSNR/sub-"$ids"/run-"$r"_tSNRmap.nii.gz $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_task-namingslowish_acq-MEMB_rec-t2star_run-"$r"_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz 

done

# calculate same statistics across runs
fslmerge -t $dirp/derivatives/tSNR/sub-"$ids"/tmp.nii.gz $dirp/derivatives/tSNR/sub-"$ids"/run-01_mean.nii.gz $dirp/derivatives/tSNR/sub-"$ids"/run-02_mean.nii.gz $dirp/derivatives/tSNR/sub-"$ids"/run-03_mean.nii.gz $dirp/derivatives/tSNR/sub-"$ids"/run-04_mean.nii.gz
3dTstat -overwrite -prefix $dirp/derivatives/tSNR/sub-"$ids"/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-"$ids"/tmp.nii.gz

fslmerge -t $dirp/derivatives/tSNR/sub-"$ids"/tmp.nii.gz $dirp/derivatives/tSNR/sub-"$ids"/run-01_stdev.nii.gz $dirp/derivatives/tSNR/sub-"$ids"/run-02_stdev.nii.gz $dirp/derivatives/tSNR/sub-"$ids"/run-03_stdev.nii.gz $dirp/derivatives/tSNR/sub-"$ids"/run-04_stdev.nii.gz
3dTstat -overwrite -prefix $dirp/derivatives/tSNR/sub-"$ids"/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-"$ids"/tmp.nii.gz

fslmerge -t $dirp/derivatives/tSNR/sub-"$ids"/tmp.nii.gz $dirp/derivatives/tSNR/sub-"$ids"/run-01_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-"$ids"/run-02_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-"$ids"/run-03_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-"$ids"/run-04_tSNRmap.nii.gz
3dTstat -overwrite -prefix $dirp/derivatives/tSNR/sub-"$ids"/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-"$ids"/tmp.nii.gz

rm $dirp/derivatives/tSNR/sub-"$ids"/tmp.nii.gz

