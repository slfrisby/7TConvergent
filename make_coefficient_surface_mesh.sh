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

# run only if volumes containing coefficients exist
if [ -d $work/sub-"$ids"/coefficients/volume/ ]; then

# make output directory
rm -rf $work/sub-"$ids"/coefficients/surface/
mkdir $work/sub-"$ids"/coefficients/surface

# get a list of the volumes
cd $work/sub-"$ids"/coefficients/volume
dir > $work/sub-"$ids"/coefficients/volume.txt

# for each volume (see "done" statement to see how input is piped in from text)
while read volume; do

# get filename for reading and saving
filename=$(basename $volume .nii)

# project to the pial surface, separately for each hemisphere
wb_command -volume-to-surface-mapping $work/sub-"$ids"/coefficients/volume/"$filename".nii $work/coefficients/pial_left.gii $work/sub-"$ids"/coefficients/surface/L_"$filename".func.gii -trilinear
wb_command -volume-to-surface-mapping $work/sub-"$ids"/coefficients/volume/"$filename".nii $work/coefficients/pial_right.gii $work/sub-"$ids"/coefficients/surface/R_"$filename".func.gii -trilinear

# smooth (6mm to match univariate - can be changed!). Treat zeros as missing data
wb_command -metric-smoothing $work/coefficients/pial_left.gii $work/sub-"$ids"/coefficients/surface/L_"$filename".func.gii 6 $work/sub-"$ids"/coefficients/surface/L_"$filename".func.gii -fix-zeros
wb_command -metric-smoothing $work/coefficients/pial_right.gii $work/sub-"$ids"/coefficients/surface/R_"$filename".func.gii 6 $work/sub-"$ids"/coefficients/surface/R_"$filename".func.gii -fix-zeros

done < $work/sub-"$ids"/coefficients/volume.txt

fi
 

