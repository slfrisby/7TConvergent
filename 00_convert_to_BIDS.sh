#!/bin/bash
export PATH=$PATH:/imaging/local/software/anaconda/latest/x86_64/bin/
# we need to access dcm2niix and it is found in
export PATH=$PATH:/group/mlr-lab/AH/Projects/AHalai
#activate FSL
FSLDIR=/imaging/local/software/fsl/latest/x86_64/fsl
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH
FSLOUTPUTTYPE=NIFTI_GZ

#conda environment with heudiconv pip installed
source activate DAISY

dirp=/imaging/projects/cbu/wbic-p00591-DAISY/main
work=/imaging/projects/cbu/wbic-p00591-DAISY/main/work
cd $dirp
# comment out the below line when using sbatch to call the script directly (rather than from within the functional preprocessing script)
ids=$1

# -p (--parents): parents or path, will also create all directories leading up to the given directory that do not exist already. For example, mkdir -p a/b will create directory a if it doesn't exist, then will create directory b inside directory a. If the given directory already exists, ignore the error.
mkdir -p $work/sub-"$ids"/

# The first time you do the analysis, you should run this line. "-c none" means that it does not carry out the conversion. This will output the tsv file dicominfo, which is found in /imaging/projects/cbu/wbic-p00591-DAISY/main/data/.heudiconv/001/info - ctrl+h to reveal hidden files on linux). Check these details against the heuristics_main.py script before running the next line. (It also produces its own heuristic.py script, which you can use as a template for the heuristics_main.py script should you need to remake or change it.)
# heudiconv -d $PWD/DICOM/sub-{subject}/*/*/*/*.dcm -o $PWD/data/ -f convertall -s 001 -c none -b --overwrite

# this line does the actual BIDS conversion. Thereafter the script denoises the MP2RAGE. It does this using the first and second inversion sequences (which are usually auto-combined by the scanner to make the MP2RAGE but luckily the originals are provided too). 
heudiconv -d $dirp/DICOM/sub-{subject}/*/*/*/*.dcm -o $dirp/data/ -f $dirp/scripts/heuristics_main.py -s "$ids" -c dcm2niix -b --overwrite

# copy unified structural files to work folder
cp $dirp/data/sub-"$ids"/anat/sub-"$ids"_run-01_T1w.nii.gz $work/sub-"$ids"/sub-"$ids"_run-01_UNI.nii.gz
cp $dirp/data/sub-"$ids"/anat/sub-"$ids"_run-01_T1w.json $work/sub-"$ids"/sub-"$ids"_run-01_UNI.json
# coregister first inversion to unified image and store this in the work folder. The transform that we apply (-applyxfm) is contained in the header of the -ref image (-usesqform expresses this.). The sort of interpolation matters only when you are changing the voxel size - nearest neighbour will do in this case.
flirt -in $dirp/data/sub-"$ids"/anat/sub-"$ids"_run-01_INV1.nii.gz  -out $work/sub-"$ids"/sub-"$ids"_run-01_INV1.nii.gz -ref $work/sub-"$ids"/sub-"$ids"_run-01_UNI.nii.gz -applyxfm -usesqform -interp nearestneighbour
# copy its .json file to the work folder
cp $dirp/data/sub-"$ids"/anat/sub-"$ids"_run-01_INV1.json $work/sub-"$ids"/sub-"$ids"_run-01_INV1.json
# do the same for the second inversion
flirt -in $dirp/data/sub-"$ids"/anat/sub-"$ids"_run-01_INV2.nii.gz  -out $work/sub-"$ids"/sub-"$ids"_run-01_INV2.nii.gz -ref $work/sub-"$ids"/sub-"$ids"_run-01_UNI.nii.gz -applyxfm -usesqform -interp nearestneighbour
cp $dirp/data/sub-"$ids"/anat/sub-"$ids"_run-01_INV2.json $work/sub-"$ids"/sub-"$ids"_run-01_INV2.json
# unzip all .nii.gz files in work folder
gzip -d $work/sub-"$ids"/*.nii.gz

#run MP2RAGE processing in MATLAB with CAT12 (using OBrian Regularisation [used] and UNI*INV2 [not used])
#requires rootdir, sub id 

matlab_r2019a -nodisplay -nodesktop -r "addpath(genpath('$dirp/scripts/'));mp2rage_proc('"$dirp"','"$ids"');exit" 

#tidy files (first and second inversions, and version unified by the scanner)
rm -rf $dirp/data/sub-"$ids"/anat/*INV* $work/sub-"$ids"/*INV* $work/sub-"$ids"/*UNI* 

# zip the outputs of CAT12
gzip $work/sub-"$ids"/mri/*.nii
#copy denoised MP2RAGE corrected image into BIDS anat folder, overwriting the version with noisy background
cp $work/sub-"$ids"/mri/msub-"$ids"_run-01_UNI_denoised.nii.gz $dirp/data/sub-"$ids"/anat/sub-"$ids"_run-01_T1w.nii.gz

#tidy fmap files - removing reverse phase files for 2nd and 3rd echoes
rm -rf $dirp/data/sub-"$ids"/fmap/*echo-2* $dirp/data/sub-"$ids"/fmap/*echo-3*

# for each run
for r in 01 02 03 04; do
# rename reverse phase files (.json and .nii.gz) to be BIDS-compliant
mv $dirp/data/sub-"$ids"/fmap/sub-"$ids"_task-namingslowish_acq-MEMB_dir-PA_run-"$r"_echo-1_epi.json $dirp/data/sub-"$ids"/fmap/sub-"$ids"_acq-MEMB_dir-PA_run-"$r"_epi.json
mv $dirp/data/sub-"$ids"/fmap/sub-"$ids"_task-namingslowish_acq-MEMB_dir-PA_run-"$r"_echo-1_epi.nii.gz $dirp/data/sub-"$ids"/fmap/sub-"$ids"_acq-MEMB_dir-PA_run-"$r"_epi.nii.gz
# get first 5 volumes (index from 0) of functional data and add that to the fmap folder
fslselectvols -i $dirp/data/sub-"$ids"/func/sub-"$ids"_task-namingslowish_acq-MEMB_run-"$r"_echo-1_bold.nii.gz -o $dirp/data/sub-"$ids"/fmap/sub-"$ids"_acq-MEMB_dir-AP_run-"$r"_epi.nii.gz --vols=0,1,2,3,4
# also copy the .json file across
cp $dirp/data/sub-"$ids"/func/sub-"$ids"_task-namingslowish_acq-MEMB_run-"$r"_echo-1_bold.json $dirp/data/sub-"$ids"/fmap/sub-"$ids"_acq-MEMB_dir-AP_run-"$r"_epi.json
done






