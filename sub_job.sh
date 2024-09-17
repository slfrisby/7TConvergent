#!/bin/bash


dirp=/imaging/projects/cbu/wbic-p00591-DAISY/main
if [ ! -d $dirp/work/logs/ ]; then
mkdir -p $dirp/work/logs/
fi

#for s in 001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024 025 026 027 028 029 030 031 032; do
for s in 023 024 025 026 027 028 029 030 031 032; do
#for s in 001; do

echo "$s"

# run preprocessing (functional script calls BIDS conversion and anatomical preprocessing if these have not already been done)
# sbatch -o $dirp/work/logs/"$s"halaiprep.out -c 16 --job-name=7T_"$s" --export=ids=${s} $dirp/scripts/02_func_preproc.sh

# runs 1st-level GLMs
#sbatch -o $dirp/work/logs/"$s"_1stglm.out -c 16 --job-name=GLM"$s" --export=ids=${s} $dirp/scripts/sub_matlabjob.sh

#transforms data into Cox format, ready for CHTC 
#sbatch -o $dirp/work/logs/"$s"_cox.out -c 16 -t 7-0:00 --job-name=COX"$s" --export=ids=${s} $dirp/scripts/sub_matlabjob.sh
#sbatch -o $dirp/work/logs/"$s"_cox.out -c 16 -q lopri -t 7-0:00 --job-name=COX"$s" --export=ids=${s} $dirp/scripts/sub_matlabjob.sh

# fits logistic LASSO models
#sbatch -o $dirp/work/logs/"$s"_LASSO.out -c 16 --job-name=LASSO"$s" --export=ids=${s} $dirp/scripts/sub_matlabjob.sh
#sbatch -o $dirp/work/logs/"$s"_LASSO.out -c 16 -q lopri --job-name=LASSO"$s" --export=ids=${s} $dirp/scripts/sub_matlabjob.sh

# RSL with LASSO
#sbatch -o $dirp/work/logs/"$s"_linear.out -c 16 --job-name=LIN"$s" --export=ids=${s} $dirp/scripts/sub_matlabjob.sh
#sbatch -o $dirp/work/logs/"$s"_linear.out -c 16 -q lopri --job-name=LIN"$s" --export=ids=${s} $dirp/scripts/sub_matlabjob.sh

# gets coefficients for LASSO models - logistic and linear
#sbatch -o $dirp/work/logs/"$s"_getcoefs.out -c 16 --job-name=COEF"$s" --export=ids=${s} $dirp/scripts/sub_matlabjob.sh
#sbatch -o $dirp/work/logs/"$s"_getcoefs.out -c 16 -q lopri --job-name=COEF"$s" --export=ids=${s} $dirp/scripts/sub_matlabjob.sh

# makes coefficient surface meshes for plotting
sbatch -o $dirp/work/logs/"$s"surface.out -c 16 --job-name=COEF_"$s" --export=ids=${s} $dirp/scripts/make_coefficient_surface_mesh.sh
#sbatch -o $dirp/work/logs/"$s"surface.out -c 16 -q lopri --job-name=COEF_"$s" --export=ids=${s} $dirp/scripts/make_coefficient_surface_mesh.sh

# models permutation distributions
# sbatch -o $dirp/work/logs/"$s"_permdist.out -c 16 --job-name=PERM"$s" --export=ids=${s} $dirp/scripts/sub_matlabjob.sh
#sbatch -o $dirp/work/logs/"$s"_permdist.out -c 16 -q lopri --job-name=PERM"$s" --export=ids=${s} $dirp/scripts/sub_matlabjob.sh

# makes tSNR maps
#sbatch -o $dirp/work/logs/"$s"tSNR.out -c 16 --job-name=7T_"$s" --export=ids=${s} $dirp/scripts/tSNR.sh

done

#001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024 025 026 027 028 029 030 031 032
	

