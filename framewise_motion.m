% by Saskia. Extracts and plots parameters necessary for making decisions
% about framewise motion.
% Based on https://github.com/FCP-INDI/C-PAC/blob/main/CPAC/generate_motion_statistics/generate_motion_statistics.py
% PARAMETERS CALCULATED:
% - Framewise displacement as per Power et al. (2012)
% - DVARS. DVARS (D temporal derivative of timecourses, VARS referring to RMS variance over voxels) indexes
%       the rate of change of BOLD signal across the entire brain at each frame of data.To calculate
%       DVARS, the volumetric timeseries is differentiated (by backwards differences) and RMS signal
%       change is calculated over the whole brain.DVARS is thus a measure of how much the intensity
%       of a brain image changes in comparison to the previous timepoint (as opposed to the global
%       signal, which is the average value of a brain image at a
%       timepoint).
% REFERENCES
% Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Spurious
%            but systematic correlations in functional connectivity MRI networks arise from subject motion. NeuroImage, 59(3),
%            2142-2154. doi:10.1016/j.neuroimage.2011.10.018

% add SPM since we need to use an HRF template
addpath('/group/mlr-lab/AH/Projects/spm12/')

% root folder
root='/imaging/projects/cbu/wbic-p00591-DAISY/main/';
cd(root);

% get all folders starting with "sub-"
folders=dir([root,'/derivatives/halaiprep/sub-*']);
% index the folders from participants with whole datasets
goodsubj = logical([1 1 1 1 0 0 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 1 1 1 1 1]);
folders = folders(goodsubj);

% initialise output. Rows are participants and columns are trials (400 for each participant with 100 for
% each run).
fd_output = zeros(size(folders,1),400);
dvars_output = zeros(size(folders,1),400);
fd_outlier = zeros(size(folders,1),1);
dvars_outlier = zeros(size(folders,1),1);
fd_thresh = zeros(size(folders,1),1);
dvars_thresh = zeros(size(folders,1),1);

for s=1:length(folders)
    % cd to folder
    cd([root,'/derivatives/halaiprep/',folders(s).name,'/func/']);

    % for each run
    for r = 1:4

        % 1. Framewise displacement as per Power et al. (2012)

        % load motion files for all 4 runs
        motion = load([folders(s).name,'_run-0',num2str(r),'_motion.txt']);
        % calculate framewise rotations and translations
        rotations = abs(diff(motion(:,1:3)));
        translations = abs(diff(motion(:,4:6)));

        % calculate FD using formula in Power paper. N.B. "Rotational
        % displacements were converted from degrees to millimeters by
        % calculating displacement on the surface of a sphere of radius
        % 50 mm, which is approximately the mean distance from the
        % cerebral cortex to the center of the head."
        fd = sum(translations, 2) + (50 * pi / 180) * sum(rotations, 2);
        % the first volume should have a FD of 0
        fd = [0; fd];

        % 2. DVARS

        % load DVARS file
        dvars = load([folders(s).name,'_run-0',num2str(r),'_DVARS.txt']);
        % the first volume has a DVARS of 0
        dvars = [0; dvars];

        % 3. Calculate motion within each trial

        % load design matrix variables
        load([root,'work/design_matrices/',folders(s).name,'_run-0',num2str(r),'_type-multivariate_design-matrix.mat']);
        
        % Actually, it isn't movement at the time given in 'onsets' that matters. It's
        % movement while we are trying to estimate the HRF. So generate an
        % HRF
        canonical_hrf = spm_hrf(1.51); % TR = 1.51
        % convolve the canonical HRF with *each trial individually*. This
        % code from http://web.mit.edu/hst.583/www/LABS/lab5b_manual.pdf
        for stim = 1:length(onsets)
            % initialise stick function
            stick_function = zeros(1,551);
            % find the index of the volume that is being acquired when that
            % stimulus is presented. Set it to 1
            trial_times = round(cell2mat(onsets(stim))/1.51);
            stick_function(trial_times)=1;
            % convolve the stick function with the HRF
            convolved_stick_function=conv(stick_function,canonical_hrf);
            convolved_stick_function= convolved_stick_function(1:551);
            % find max value of HRF
            [m, idx] = max(convolved_stick_function);
            % get fd data for a window spanning 2 volumes either side of
            % the max(5 vols in total)
            fd_window = fd(idx-2:idx+2);
            % average
            fd_mean = mean(fd_window);
            % store
            fd_output(s,(100*(r-1)+stim)) = fd_mean;
            % do the same for dvars
            dvars_window = dvars(idx-2:idx+2);
            dvars_mean = mean(dvars_window);
            dvars_output(s,(100*(r-1)+stim)) = dvars_mean;
        end
    end
    % calculate 2 metrics of quality:
    % 1. the number of trials with FD or DVARS greater than 2 standard
    % deviations above the mean
    % 2. the number of trials with FD or DVARS above the thresholds
    % implemented in fmriprep (0.5 for FD, 1.5 for DVARS)
    fd_outlier(s) = length(find(fd_output(s,:)>mean(fd_output(s,:))+2*std(fd_output(s,:))));
    dvars_outlier(s) = length(find(dvars_output(s,:)>mean(dvars_output(s,:))+2*std(dvars_output(s,:))));
    fd_thresh(s) = length(find(fd_output(s,:)>0.5));
    dvars_thresh(s) = length(find(dvars_output(s,:)>1.5));
end

figure;
% plot absolute FD and DVARS
subplot(3,2,1)
boxplot(fd_output')
title('Framewise Displacement')
xticklabels({folders.name})
% plot line for fmriprep cutoff
yline(0.5,'r:')
subplot(3,2,2)
boxplot(dvars_output')
title('DVARS')
xticklabels({folders.name})
yline(1.5,'r:')
% plot number of trials above threshold
subplot(3,2,3)
bar({folders.name},fd_thresh','r')
title('Framewise Displacement - n above threshold')
subplot(3,2,4)
bar({folders.name},dvars_thresh','b')
title('DVARS - n above threshold')
% plot number of outliers
subplot(3,2,5)
bar({folders.name},fd_outlier','r')
title('Framewise Displacement - outliers')
subplot(3,2,6)
bar({folders.name},dvars_outlier','b')
title('DVARS - outliers')