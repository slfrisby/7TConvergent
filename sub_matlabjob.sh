#!/bin/bash
echo "
++++++++++++++++++++++++" 

#ids is set in sub_job.sh 
#ids=$1
dirp=/imaging/projects/cbu/wbic-p00591-DAISY/main
work=/imaging/projects/cbu/wbic-p00591-DAISY/main/work

# runs 1st level GLMs
#matlab_r2023b -nodisplay -nodesktop -r "addpath('$dirp/scripts/');firstlevel_glm('"$ids"');exit"

# does Cox formatting
#matlab_r2023b -nodisplay -nodesktop -r "addpath('$dirp/scripts/');cox_format('"$ids"');exit"

# fits LASSO models for classification
#matlab_r2023b -nodisplay -nodesktop -r "addpath('$dirp/scripts/');fit_log_models('"$ids"');exit"

# fits linear models with LASSO to predict each dimension
matlab_r2023b -nodisplay -nodesktop -r "addpath('$dirp/scripts/');fit_linear_models('"$ids"');exit"

# gets coefficients for LASSO models - logistic and linear
#matlab_r2023b -nodisplay -nodesktop -r "addpath('$dirp/scripts/');get_LASSO_coefficients('"$ids"');exit"

# models permutation distributions
# matlab_r2023b -nodisplay -nodesktop -r "addpath('$dirp/scripts/');generate_permutation_distribution('"$ids"');exit"

