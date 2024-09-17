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
cd $dirp

# concatenate maps from each participant (except excluded particpants) - mean, standard deviation, tSNR - then calculate mean over participants
fslmerge -t $work/group_mean.nii.gz $dirp/derivatives/tSNR/sub-001/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-002/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-003/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-004/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-007/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-009/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-010/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-011/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-012/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-013/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-014/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-015/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-016/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-017/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-018/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-019/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-020/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-021/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-022/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-023/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-024/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-026/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-028/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-029/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-030/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-031/acrossrun_mean.nii.gz $dirp/derivatives/tSNR/sub-032/acrossrun_mean.nii.gz

3dTstat -overwrite -prefix $dirp/derivatives/tSNR/group_mean.nii.gz $work/group_mean.nii.gz

fslmerge -t $work/group_stdev.nii.gz $dirp/derivatives/tSNR/sub-001/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-002/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-003/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-004/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-007/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-009/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-010/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-011/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-012/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-013/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-014/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-015/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-016/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-017/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-018/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-019/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-020/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-021/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-022/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-023/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-024/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-026/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-028/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-029/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-030/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-031/acrossrun_stdev.nii.gz $dirp/derivatives/tSNR/sub-032/acrossrun_stdev.nii.gz

3dTstat -overwrite -prefix $dirp/derivatives/tSNR/group_stdev.nii.gz $work/group_stdev.nii.gz

fslmerge -t $work/group_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-001/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-002/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-003/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-004/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-007/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-009/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-010/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-011/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-012/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-013/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-014/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-015/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-016/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-017/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-018/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-019/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-020/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-021/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-022/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-023/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-024/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-026/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-028/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-029/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-030/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-031/acrossrun_tSNRmap.nii.gz $dirp/derivatives/tSNR/sub-032/acrossrun_tSNRmap.nii.gz

3dTstat -overwrite -prefix $dirp/derivatives/tSNR/group_tSNRmap.nii.gz $work/group_tSNRmap.nii.gz

