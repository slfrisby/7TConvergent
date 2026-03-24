% collates results and plots figures (an updated version of
% paper_bar_charts.m). 

% setup
addpath('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/');
addpath('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/WISC_MVPA/src/');
root = '/imaging/projects/cbu/wbic-p00591-DAISY/main/derivatives';
cd(root);

% load metadata for all participants - this provides cross-validation
% indices
load([root,'/cox/tedana_averaged_metadata.mat']);

% specify IDs of participants to include in analysis
subcode = {'sub-001','sub-002','sub-003','sub-004','sub-007','sub-009','sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019','sub-020','sub-021','sub-022','sub-023','sub-024','sub-026','sub-028','sub-029','sub-030','sub-031','sub-032'};
% specify ROIs
ROI = {'wholebrain','LATL','LATLant','LATLpos','RATL','RATLant','RATLpos'};

% get target dimensions for RSL
load('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/dilkina_norms.mat');
[C,z] = embed_similarity_matrix(dilkina_norms,3);
U = rescale_embedding(C,z);

% define colour palettes
colours = load(['/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/colour_palettes.mat']);

% organise results from logistic regression with LASSO. These were
% generated with handwritten code run on the CBU cluster.

% for every participant
for s = 1:size(subcode,2)

    % load results
    tmp = load([root,'/LASSO/log/tedana/',subcode{s},'/results.mat']);
    
    % organise main results - fold x participant
    logLASSO.wholebrain.final(:,s) = tmp.averaged_whole_brain.output;
    logLASSO.LATL.final(:,s) = tmp.averaged_left_ATL.output;
    logLASSO.LATLant.final(:,s) = tmp.averaged_left_ATL_ant.output;
    logLASSO.LATLpos.final(:,s) = tmp.averaged_left_ATL_pos.output;
    logLASSO.RATL.final(:,s) = tmp.averaged_right_ATL.output;
    logLASSO.RATLant.final(:,s) = tmp.averaged_right_ATL_ant.output;
    logLASSO.RATLpos.final(:,s) = tmp.averaged_right_ATL_pos.output;

    % organise permutation results: seed x participant
    logLASSO.wholebrain.perm(:,s) = mean(tmp.averaged_whole_brain.permoutput)';
    logLASSO.LATL.perm(:,s) = mean(tmp.averaged_left_ATL.permoutput)';
    logLASSO.LATLant.perm(:,s) = mean(tmp.averaged_left_ATL_ant.permoutput)';
    logLASSO.LATLpos.perm(:,s) = mean(tmp.averaged_left_ATL_pos.permoutput)';
    logLASSO.RATL.perm(:,s) = mean(tmp.averaged_right_ATL.permoutput)';
    logLASSO.RATLant.perm(:,s) = mean(tmp.averaged_right_ATL_ant.permoutput)';
    logLASSO.RATLpos.perm(:,s) = mean(tmp.averaged_right_ATL_pos.permoutput)';
end

% then calculate the mean over holdout folds
for r = 1:size(ROI,2)
    % calculate mean over holdout folds
    logLASSO.(ROI{r}).final = mean(logLASSO.(ROI{r}).final,1);
end

% organise results from logistic regression with SOSLASSO. These were
% generated with the WISC_MVPA toolbox run on CHTC in Madison.

% for each ROI
for r = 1:size(ROI,2)

    % load main results
    load([root,'/condor/derivatives/classification/SOSLASSO/performance/',ROI{r},'/final_performance.mat']);

    % for each participant
    for s = 1:size(subcode,2)
        subjectIndex = str2num(erase(subcode{s},'sub-'));
        % for each holdout fold
        for ho = 1:10
            % get results
            tmp = Tallcv((Tallcv.subject == subjectIndex & Tallcv.cvholdout == ho), :);
            % get predicted coordinates for that holdout fold
            tmp = tmp.Yz{1};
            predictedCoords = tmp(metadata(s).cvind(:,1) == ho);
            % compare with target (0 for animate, 1 for inanimate) and get
            % AUC
            [~, ~, ~, AUC] = perfcurve([zeros(5,1);ones(5,1)], predictedCoords, 1);
            % organise
            SOSLASSO.(ROI{r}).final(ho,s) = AUC;
        end    
    end

    % load permutation results
    load([root,'/condor/derivatives/classification/SOSLASSO/performance/',ROI{r},'/perm_performance.mat']);

    % for each participant
    for s = 1:size(subcode,2)
        subjectIndex = str2num(erase(subcode{s},'sub-'));
        % for each random seed
        for seed = 1:100 
            % for each holdout set
            for ho = 1:10
                % get results
                tmp = Tallcv((Tallcv.subject == subjectIndex & Tallcv.cvholdout == ho & Tallcv.RandomSeed == seed), :);
                % get predicted coordinates for that holdout fold
                tmp = tmp.Yz{1};
                predictedCoords = tmp(metadata(s).cvind(:,1) == ho);
                % compare with target (0 for animate, 1 for inanimate) and get
                % AUC
                [~, ~, ~, AUC] = perfcurve([zeros(5,1);ones(5,1)], predictedCoords, 1);
                SOSLASSO.(ROI{r}).perm(seed,s,ho) = AUC;
            end
        end
    end
    % calculate mean over holdout folds
    SOSLASSO.(ROI{r}).final = mean(SOSLASSO.(ROI{r}).final,1);
    SOSLASSO.(ROI{r}).perm = mean(SOSLASSO.(ROI{r}).perm,3);
end

% organise results from RSL with LASSO. These were
% generated with handwritten code run on the CBU cluster.

for s = 1:size(subcode,2)

    % load results
    tmp = load([root,'/LASSO/linear/tedana/',subcode{s},'/results.mat']);
    
    % organise main results: fold x participant
    for ho = 1:10
        % get target coordinates for D1
        targetCoords = U(metadata(s).cvind(:,1) == ho, 1);
        % get predicted coordinates
        predictedCoords = tmp.averaged_whole_brain_D1.predictedcoords(metadata(s).cvind(:,1) == ho);
        % calculate correlation across all stimuli
        linearLASSO.wholebrain.D1.all.final(ho,s) = corr(targetCoords,predictedCoords);
        % calculate correlation for just animate and just inanimate stimuli
        linearLASSO.wholebrain.D1.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.wholebrain.D1.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        % repeat for all ROIs
        predictedCoords = tmp.averaged_left_ATL_D1.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.LATL.D1.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.LATL.D1.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.LATL.D1.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_left_ATL_ant_D1.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.LATLant.D1.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.LATLant.D1.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.LATLant.D1.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_left_ATL_pos_D1.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.LATLpos.D1.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.LATLpos.D1.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.LATLpos.D1.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_right_ATL_D1.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.RATL.D1.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.RATL.D1.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.RATL.D1.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_right_ATL_ant_D1.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.RATLant.D1.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.RATLant.D1.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.RATLant.D1.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_right_ATL_pos_D1.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.RATLpos.D1.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.RATLpos.D1.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.RATLpos.D1.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        % and for dimension 2
        targetCoords = U(metadata(s).cvind(:,1) == ho, 2);
        predictedCoords = tmp.averaged_whole_brain_D2.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.wholebrain.D2.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.wholebrain.D2.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.wholebrain.D2.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_left_ATL_D2.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.LATL.D2.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.LATL.D2.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.LATL.D2.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_left_ATL_ant_D2.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.LATLant.D2.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.LATLant.D2.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.LATLant.D2.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_left_ATL_pos_D2.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.LATLpos.D2.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.LATLpos.D2.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.LATLpos.D2.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_right_ATL_D2.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.RATL.D2.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.RATL.D2.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.RATL.D2.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_right_ATL_ant_D2.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.RATLant.D2.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.RATLant.D2.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.RATLant.D2.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_right_ATL_pos_D2.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.RATLpos.D2.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.RATLpos.D2.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.RATLpos.D2.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        % and for dimension 3
        targetCoords = U(metadata(s).cvind(:,1) == ho, 3);
        predictedCoords = tmp.averaged_whole_brain_D3.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.wholebrain.D3.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.wholebrain.D3.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.wholebrain.D3.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_left_ATL_D3.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.LATL.D3.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.LATL.D3.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.LATL.D3.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_left_ATL_ant_D3.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.LATLant.D3.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.LATLant.D3.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.LATLant.D3.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_left_ATL_pos_D3.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.LATLpos.D3.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.LATLpos.D3.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.LATLpos.D3.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_right_ATL_D3.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.RATL.D3.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.RATL.D3.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.RATL.D3.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_right_ATL_ant_D3.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.RATLant.D3.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.RATLant.D3.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.RATLant.D3.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
        predictedCoords = tmp.averaged_right_ATL_pos_D3.predictedcoords(metadata(s).cvind(:,1) == ho);
        linearLASSO.RATLpos.D3.all.final(ho,s) = corr(targetCoords,predictedCoords);
        linearLASSO.RATLpos.D3.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
        linearLASSO.RATLpos.D3.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
    end
     
    % organise permutation results: seed x participant
    for seed = 1:100
        % initialise output
        for ho = 1:10
            % get target coordinates for D1
            targetCoords = U(metadata(s).cvind(:,1) == ho, 1);
            % get predicted coordinates
            predictedCoords = tmp.averaged_whole_brain_D1.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            % calculate correlation across all stimuli
            linearLASSO.wholebrain.D1.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            % calculate correlation for just animate and just inanimate stimuli
            linearLASSO.wholebrain.D1.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.wholebrain.D1.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            % repeat for all ROIs
            predictedCoords = tmp.averaged_left_ATL_D1.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.LATL.D1.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.LATL.D1.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.LATL.D1.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_left_ATL_ant_D1.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.LATLant.D1.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.LATLant.D1.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.LATLant.D1.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_left_ATL_pos_D1.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.LATLpos.D1.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.LATLpos.D1.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.LATLpos.D1.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_right_ATL_D1.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.RATL.D1.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.RATL.D1.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.RATL.D1.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_right_ATL_ant_D1.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.RATLant.D1.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.RATLant.D1.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.RATLant.D1.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_right_ATL_pos_D1.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.RATLpos.D1.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.RATLpos.D1.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.RATLpos.D1.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            % and for dimension 2
            targetCoords = U(metadata(s).cvind(:,1) == ho, 2);
            predictedCoords = tmp.averaged_whole_brain_D2.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.wholebrain.D2.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.wholebrain.D2.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.wholebrain.D2.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_left_ATL_D2.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.LATL.D2.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.LATL.D2.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.LATL.D2.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_left_ATL_ant_D2.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.LATLant.D2.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.LATLant.D2.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.LATLant.D2.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_left_ATL_pos_D2.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.LATLpos.D2.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.LATLpos.D2.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.LATLpos.D2.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_right_ATL_D2.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.RATL.D2.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.RATL.D2.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.RATL.D2.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_right_ATL_ant_D2.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.RATLant.D2.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.RATLant.D2.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.RATLant.D2.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_right_ATL_pos_D2.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.RATLpos.D2.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.RATLpos.D2.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.RATLpos.D2.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            % and for dimension 3
            targetCoords = U(metadata(s).cvind(:,1) == ho, 3);
            predictedCoords = tmp.averaged_whole_brain_D3.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.wholebrain.D3.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.wholebrain.D3.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.wholebrain.D3.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_left_ATL_D3.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.LATL.D3.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.LATL.D3.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.LATL.D3.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_left_ATL_ant_D3.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.LATLant.D3.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.LATLant.D3.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.LATLant.D3.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_left_ATL_pos_D3.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.LATLpos.D3.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.LATLpos.D3.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.LATLpos.D3.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_right_ATL_D3.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.RATL.D3.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.RATL.D3.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.RATL.D3.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_right_ATL_ant_D3.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.RATLant.D3.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.RATLant.D3.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.RATLant.D3.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
            predictedCoords = tmp.averaged_right_ATL_pos_D3.permpredictedcoords(metadata(s).cvind(:,1) == ho, seed);
            linearLASSO.RATLpos.D3.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
            linearLASSO.RATLpos.D3.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
            linearLASSO.RATLpos.D3.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
        end
    end
end

% in some cases the model predicts the same value for everything. If this
% is the case, the correlation is zero but (since the variance is zero) it
% will be recorded as NaN. Correct this using a recursive function
linearLASSO = replace_nan(linearLASSO);

% then calculate the mean over holdout folds
for r = 1:size(ROI,2)
    % calculate mean over holdout folds
    linearLASSO.(ROI{r}).D1.all.final = mean(linearLASSO.(ROI{r}).D1.all.final,1);
    linearLASSO.(ROI{r}).D1.animate.final = mean(linearLASSO.(ROI{r}).D1.animate.final,1);
    linearLASSO.(ROI{r}).D1.inanimate.final = mean(linearLASSO.(ROI{r}).D1.inanimate.final,1);
    linearLASSO.(ROI{r}).D2.all.final = mean(linearLASSO.(ROI{r}).D2.all.final,1);
    linearLASSO.(ROI{r}).D2.animate.final = mean(linearLASSO.(ROI{r}).D2.animate.final,1);
    linearLASSO.(ROI{r}).D2.inanimate.final = mean(linearLASSO.(ROI{r}).D2.inanimate.final,1);
    linearLASSO.(ROI{r}).D3.all.final = mean(linearLASSO.(ROI{r}).D3.all.final,1);
    linearLASSO.(ROI{r}).D3.animate.final = mean(linearLASSO.(ROI{r}).D3.animate.final,1);
    linearLASSO.(ROI{r}).D3.inanimate.final = mean(linearLASSO.(ROI{r}).D3.inanimate.final,1);
    linearLASSO.(ROI{r}).D1.all.perm = mean(linearLASSO.(ROI{r}).D1.all.perm,3);
    linearLASSO.(ROI{r}).D1.animate.perm = mean(linearLASSO.(ROI{r}).D1.animate.perm,3);
    linearLASSO.(ROI{r}).D1.inanimate.perm = mean(linearLASSO.(ROI{r}).D1.inanimate.perm,3);
    linearLASSO.(ROI{r}).D2.all.perm = mean(linearLASSO.(ROI{r}).D2.all.perm,3);
    linearLASSO.(ROI{r}).D2.animate.perm = mean(linearLASSO.(ROI{r}).D2.animate.perm,3);
    linearLASSO.(ROI{r}).D2.inanimate.perm = mean(linearLASSO.(ROI{r}).D2.inanimate.perm,3);
    linearLASSO.(ROI{r}).D3.all.perm = mean(linearLASSO.(ROI{r}).D3.all.perm,3);
    linearLASSO.(ROI{r}).D3.animate.perm = mean(linearLASSO.(ROI{r}).D3.animate.perm,3);
    linearLASSO.(ROI{r}).D3.inanimate.perm = mean(linearLASSO.(ROI{r}).D3.inanimate.perm,3);
end

% organise results from RSL with grOWL. These were
% generated with the WISC_MVPA toolbox run on CHTC in Madison.

% for each ROI
for r = 1:size(ROI,2)

    % load main results
    load([root,'/condor/derivatives/correlation/grOWL/performance/',ROI{r},'/final_performance.mat']);

    % for each participant
    for s = 1:size(subcode,2)
        subjectIndex = str2num(erase(subcode{s},'sub-'));
        % for each holdout fold
        for ho = 1:10
            % get results
            tmp = Tallcv((Tallcv.subject == subjectIndex & Tallcv.cvholdout == ho), :);
            tmp = tmp.Cz{1};
            % get predicted coordinates for that holdout fold for D1
            predictedCoords = tmp(metadata(s).cvind(:,1) == ho, 1);
            % get target coordinates for D1
            targetCoords = U(metadata(s).cvind(:,1) == ho, 1);
            % calculate correlation for all stimuli
            grOWL.(ROI{r}).D1.all.final(ho,s) = corr(targetCoords,predictedCoords);
            % calculate correlation just for animate and just for inanimate stimuli
            grOWL.(ROI{r}).D1.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
            grOWL.(ROI{r}).D1.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
            % repeat for dimension 2
            predictedCoords = tmp(metadata(s).cvind(:,1) == ho, 2);
            targetCoords = U(metadata(s).cvind(:,1) == ho, 2);
            grOWL.(ROI{r}).D2.all.final(ho,s) = corr(targetCoords,predictedCoords);
            grOWL.(ROI{r}).D2.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
            grOWL.(ROI{r}).D2.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));
            % repeat for dimension 3
            predictedCoords = tmp(metadata(s).cvind(:,1) == ho, 3);
            targetCoords = U(metadata(s).cvind(:,1) == ho, 3);
            grOWL.(ROI{r}).D3.all.final(ho,s) = corr(targetCoords,predictedCoords);
            grOWL.(ROI{r}).D3.animate.final(ho,s) = corr(targetCoords(1:5),predictedCoords(1:5));
            grOWL.(ROI{r}).D3.inanimate.final(ho,s) = corr(targetCoords(6:10),predictedCoords(6:10));          
        end    
    end

    % load permutation results
    load([root,'/condor/derivatives/correlation/grOWL/performance/',ROI{r},'/perm_performance.mat']);

    % for each participant
    for s = 1:size(subcode,2)
        subjectIndex = str2num(erase(subcode{s},'sub-'));
        % for each random seed
        for seed = 1:100 
            % for each holdout set
            for ho = 1:10
                % get results
                tmp = Tallcv((Tallcv.subject == subjectIndex & Tallcv.cvholdout == ho & Tallcv.RandomSeed == seed), :);
                tmp = tmp.Cz{1};
                % get predicted coordinates for that holdout fold for D1
                predictedCoords = tmp(metadata(s).cvind(:,1) == ho, 1);
                % get target coordinates for D1
                targetCoords = U(metadata(s).cvind(:,1) == ho, 1);
                % calculate correlation for all stimuli
                grOWL.(ROI{r}).D1.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
                % calculate correlation just for animate and just for inanimate stimuli
                grOWL.(ROI{r}).D1.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
                grOWL.(ROI{r}).D1.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
                % repeat for dimension 2
                predictedCoords = tmp(metadata(s).cvind(:,1) == ho, 2);
                targetCoords = U(metadata(s).cvind(:,1) == ho, 2);
                grOWL.(ROI{r}).D2.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
                grOWL.(ROI{r}).D2.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
                grOWL.(ROI{r}).D2.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));
                % repeat for dimension 3
                predictedCoords = tmp(metadata(s).cvind(:,1) == ho, 3);
                targetCoords = U(metadata(s).cvind(:,1) == ho, 3);
                grOWL.(ROI{r}).D3.all.perm(seed,s,ho) = corr(targetCoords,predictedCoords);
                grOWL.(ROI{r}).D3.animate.perm(seed,s,ho) = corr(targetCoords(1:5),predictedCoords(1:5));
                grOWL.(ROI{r}).D3.inanimate.perm(seed,s,ho) = corr(targetCoords(6:10),predictedCoords(6:10));    
            end
        end
    end
end

% deal with nans
grOWL = replace_nan(grOWL);

% then calculate the mean over holdout folds
for r = 1:size(ROI,2)
    % calculate mean over holdout folds
    grOWL.(ROI{r}).D1.all.final = mean(grOWL.(ROI{r}).D1.all.final,1);
    grOWL.(ROI{r}).D1.animate.final = mean(grOWL.(ROI{r}).D1.animate.final,1);
    grOWL.(ROI{r}).D1.inanimate.final = mean(grOWL.(ROI{r}).D1.inanimate.final,1);
    grOWL.(ROI{r}).D2.all.final = mean(grOWL.(ROI{r}).D2.all.final,1);
    grOWL.(ROI{r}).D2.animate.final = mean(grOWL.(ROI{r}).D2.animate.final,1);
    grOWL.(ROI{r}).D2.inanimate.final = mean(grOWL.(ROI{r}).D2.inanimate.final,1);
    grOWL.(ROI{r}).D3.all.final = mean(grOWL.(ROI{r}).D3.all.final,1);
    grOWL.(ROI{r}).D3.animate.final = mean(grOWL.(ROI{r}).D3.animate.final,1);
    grOWL.(ROI{r}).D3.inanimate.final = mean(grOWL.(ROI{r}).D3.inanimate.final,1);

    grOWL.(ROI{r}).D1.all.perm = mean(grOWL.(ROI{r}).D1.all.perm,3);
    grOWL.(ROI{r}).D1.animate.perm = mean(grOWL.(ROI{r}).D1.animate.perm,3);
    grOWL.(ROI{r}).D1.inanimate.perm = mean(grOWL.(ROI{r}).D1.inanimate.perm,3);
    grOWL.(ROI{r}).D2.all.perm = mean(grOWL.(ROI{r}).D2.all.perm,3);
    grOWL.(ROI{r}).D2.animate.perm = mean(grOWL.(ROI{r}).D2.animate.perm,3);
    grOWL.(ROI{r}).D2.inanimate.perm = mean(grOWL.(ROI{r}).D2.inanimate.perm,3);
    grOWL.(ROI{r}).D3.all.perm = mean(grOWL.(ROI{r}).D3.all.perm,3);
    grOWL.(ROI{r}).D3.animate.perm = mean(grOWL.(ROI{r}).D3.animate.perm,3);
    grOWL.(ROI{r}).D3.inanimate.perm = mean(grOWL.(ROI{r}).D3.inanimate.perm,3);
end

% also generate variables containing the difference between classification with logistic LASSO
% and SOSLASSO, and between RSL with linear LASSO and with grOWL.
for r = 1:size(ROI,2)

    classificationDifference.(ROI{r}).final = SOSLASSO.(ROI{r}).final - logLASSO.(ROI{r}).final;
    classificationDifference.(ROI{r}).perm = SOSLASSO.(ROI{r}).perm - logLASSO.(ROI{r}).perm;

    RSLDifference.(ROI{r}).D1.all.final = grOWL.(ROI{r}).D1.all.final - linearLASSO.(ROI{r}).D1.all.final;
    RSLDifference.(ROI{r}).D1.animate.final = grOWL.(ROI{r}).D1.animate.final - linearLASSO.(ROI{r}).D1.animate.final;
    RSLDifference.(ROI{r}).D1.inanimate.final = grOWL.(ROI{r}).D1.inanimate.final - linearLASSO.(ROI{r}).D1.inanimate.final;
    RSLDifference.(ROI{r}).D2.all.final = grOWL.(ROI{r}).D2.all.final - linearLASSO.(ROI{r}).D2.all.final;
    RSLDifference.(ROI{r}).D2.animate.final = grOWL.(ROI{r}).D2.animate.final - linearLASSO.(ROI{r}).D2.animate.final;
    RSLDifference.(ROI{r}).D2.inanimate.final = grOWL.(ROI{r}).D2.inanimate.final - linearLASSO.(ROI{r}).D2.inanimate.final;
    RSLDifference.(ROI{r}).D3.all.final = grOWL.(ROI{r}).D3.all.final - linearLASSO.(ROI{r}).D3.all.final;
    RSLDifference.(ROI{r}).D3.animate.final = grOWL.(ROI{r}).D3.animate.final - linearLASSO.(ROI{r}).D3.animate.final;
    RSLDifference.(ROI{r}).D3.inanimate.final = grOWL.(ROI{r}).D3.inanimate.final - linearLASSO.(ROI{r}).D3.inanimate.final;
    RSLDifference.(ROI{r}).D1.all.perm = grOWL.(ROI{r}).D1.all.perm - linearLASSO.(ROI{r}).D1.all.perm;
    RSLDifference.(ROI{r}).D1.animate.perm = grOWL.(ROI{r}).D1.animate.perm - linearLASSO.(ROI{r}).D1.animate.perm;
    RSLDifference.(ROI{r}).D1.inanimate.perm = grOWL.(ROI{r}).D1.inanimate.perm - linearLASSO.(ROI{r}).D1.inanimate.perm;
    RSLDifference.(ROI{r}).D2.all.perm = grOWL.(ROI{r}).D2.all.perm - linearLASSO.(ROI{r}).D2.all.perm;
    RSLDifference.(ROI{r}).D2.animate.perm = grOWL.(ROI{r}).D2.animate.perm - linearLASSO.(ROI{r}).D2.animate.perm;
    RSLDifference.(ROI{r}).D2.inanimate.perm = grOWL.(ROI{r}).D2.inanimate.perm - linearLASSO.(ROI{r}).D2.inanimate.perm;
    RSLDifference.(ROI{r}).D3.all.perm = grOWL.(ROI{r}).D3.all.perm - linearLASSO.(ROI{r}).D3.all.perm;
    RSLDifference.(ROI{r}).D3.animate.perm = grOWL.(ROI{r}).D3.animate.perm - linearLASSO.(ROI{r}).D3.animate.perm;
    RSLDifference.(ROI{r}).D3.inanimate.perm = grOWL.(ROI{r}).D3.inanimate.perm - linearLASSO.(ROI{r}).D3.inanimate.perm;   

end

% and a variable for testing whether the difference between anterior and
% posterior ROIs is significant. 
dynamismDifference.LASSO.left.final = logLASSO.LATLpos.final - logLASSO.LATLant.final;
dynamismDifference.LASSO.left.perm = logLASSO.LATLpos.perm - logLASSO.LATLant.perm;
dynamismDifference.LASSO.right.final = logLASSO.RATLpos.final - logLASSO.RATLant.final;
dynamismDifference.LASSO.right.perm = logLASSO.RATLpos.perm - logLASSO.RATLant.perm;
dynamismDifference.SOSLASSO.left.final = SOSLASSO.LATLpos.final - SOSLASSO.LATLant.final;
dynamismDifference.SOSLASSO.left.perm = SOSLASSO.LATLpos.perm - SOSLASSO.LATLant.perm;
dynamismDifference.SOSLASSO.right.final = SOSLASSO.RATLpos.final - SOSLASSO.RATLant.final;
dynamismDifference.SOSLASSO.right.perm = SOSLASSO.RATLpos.perm - SOSLASSO.RATLant.perm;

% next, Stelzer-bootstrap permutation distributions.
for r = 1:size(ROI,2)

    logLASSO.(ROI{r}).bootstrappedPerm = stelzer_bootstrap(logLASSO.(ROI{r}).perm);
    SOSLASSO.(ROI{r}).bootstrappedPerm = stelzer_bootstrap(SOSLASSO.(ROI{r}).perm);

    classificationDifference.(ROI{r}).bootstrappedPerm = stelzer_bootstrap(classificationDifference.(ROI{r}).perm);

    linearLASSO.(ROI{r}).D1.all.bootstrappedPerm = stelzer_bootstrap(linearLASSO.(ROI{r}).D1.all.perm);
    linearLASSO.(ROI{r}).D1.animate.bootstrappedPerm = stelzer_bootstrap(linearLASSO.(ROI{r}).D1.animate.perm);
    linearLASSO.(ROI{r}).D1.inanimate.bootstrappedPerm = stelzer_bootstrap(linearLASSO.(ROI{r}).D1.inanimate.perm);
    linearLASSO.(ROI{r}).D2.all.bootstrappedPerm = stelzer_bootstrap(linearLASSO.(ROI{r}).D2.all.perm);
    linearLASSO.(ROI{r}).D2.animate.bootstrappedPerm = stelzer_bootstrap(linearLASSO.(ROI{r}).D2.animate.perm);
    linearLASSO.(ROI{r}).D2.inanimate.bootstrappedPerm = stelzer_bootstrap(linearLASSO.(ROI{r}).D2.inanimate.perm);
    linearLASSO.(ROI{r}).D3.all.bootstrappedPerm = stelzer_bootstrap(linearLASSO.(ROI{r}).D3.all.perm);
    linearLASSO.(ROI{r}).D3.animate.bootstrappedPerm = stelzer_bootstrap(linearLASSO.(ROI{r}).D3.animate.perm);
    linearLASSO.(ROI{r}).D3.inanimate.bootstrappedPerm = stelzer_bootstrap(linearLASSO.(ROI{r}).D3.inanimate.perm);
    grOWL.(ROI{r}).D1.all.bootstrappedPerm = stelzer_bootstrap(grOWL.(ROI{r}).D1.all.perm);
    grOWL.(ROI{r}).D1.animate.bootstrappedPerm = stelzer_bootstrap(grOWL.(ROI{r}).D1.animate.perm);
    grOWL.(ROI{r}).D1.inanimate.bootstrappedPerm = stelzer_bootstrap(grOWL.(ROI{r}).D1.inanimate.perm);
    grOWL.(ROI{r}).D2.all.bootstrappedPerm = stelzer_bootstrap(grOWL.(ROI{r}).D2.all.perm);
    grOWL.(ROI{r}).D2.animate.bootstrappedPerm = stelzer_bootstrap(grOWL.(ROI{r}).D2.animate.perm);
    grOWL.(ROI{r}).D2.inanimate.bootstrappedPerm = stelzer_bootstrap(grOWL.(ROI{r}).D2.inanimate.perm);
    grOWL.(ROI{r}).D3.all.bootstrappedPerm = stelzer_bootstrap(grOWL.(ROI{r}).D3.all.perm);
    grOWL.(ROI{r}).D3.animate.bootstrappedPerm = stelzer_bootstrap(grOWL.(ROI{r}).D3.animate.perm);
    grOWL.(ROI{r}).D3.inanimate.bootstrappedPerm = stelzer_bootstrap(grOWL.(ROI{r}).D3.inanimate.perm);

    RSLDifference.(ROI{r}).D1.all.bootstrappedPerm = stelzer_bootstrap(RSLDifference.(ROI{r}).D1.all.perm);
    RSLDifference.(ROI{r}).D1.animate.bootstrappedPerm = stelzer_bootstrap(RSLDifference.(ROI{r}).D1.animate.perm);
    RSLDifference.(ROI{r}).D1.inanimate.bootstrappedPerm = stelzer_bootstrap(RSLDifference.(ROI{r}).D1.inanimate.perm);
    RSLDifference.(ROI{r}).D2.all.bootstrappedPerm = stelzer_bootstrap(RSLDifference.(ROI{r}).D2.all.perm);
    RSLDifference.(ROI{r}).D2.animate.bootstrappedPerm = stelzer_bootstrap(RSLDifference.(ROI{r}).D2.animate.perm);
    RSLDifference.(ROI{r}).D2.inanimate.bootstrappedPerm = stelzer_bootstrap(RSLDifference.(ROI{r}).D2.inanimate.perm);
    RSLDifference.(ROI{r}).D3.all.bootstrappedPerm = stelzer_bootstrap(RSLDifference.(ROI{r}).D3.all.perm);
    RSLDifference.(ROI{r}).D3.animate.bootstrappedPerm = stelzer_bootstrap(RSLDifference.(ROI{r}).D3.animate.perm);
    RSLDifference.(ROI{r}).D3.inanimate.bootstrappedPerm = stelzer_bootstrap(RSLDifference.(ROI{r}).D3.inanimate.perm);
end

dynamismDifference.LASSO.left.bootstrappedPerm = stelzer_bootstrap(dynamismDifference.LASSO.left.perm);
dynamismDifference.LASSO.right.bootstrappedPerm = stelzer_bootstrap(dynamismDifference.LASSO.right.perm);
dynamismDifference.SOSLASSO.left.bootstrappedPerm = stelzer_bootstrap(dynamismDifference.SOSLASSO.left.perm);
dynamismDifference.SOSLASSO.right.bootstrappedPerm = stelzer_bootstrap(dynamismDifference.SOSLASSO.right.perm);

% make folder in which to store plots
if ~exist('/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts')
    mkdir('/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts')
end

% plot classifications. LASSO and SOSLASSO for the ATL

clear plot
% prepare the variables - mean over participants
plot.mean = [mean(logLASSO.LATL.final), mean(logLASSO.LATLant.final), mean(logLASSO.LATLpos.final); ...
    mean(logLASSO.RATL.final), mean(logLASSO.RATLant.final), mean(logLASSO.RATLpos.final); ...
    mean(SOSLASSO.LATL.final), mean(SOSLASSO.LATLant.final), mean(SOSLASSO.LATLpos.final); ...
    mean(SOSLASSO.RATL.final), mean(SOSLASSO.RATLant.final), mean(SOSLASSO.RATLpos.final)];
% standard error
plot.standardError = [(std(logLASSO.LATL.final)/sqrt(size(logLASSO.LATL.final,2))), (std(logLASSO.LATLant.final)/sqrt(size(logLASSO.LATLant.final,2))), (std(logLASSO.LATLpos.final)/sqrt(size(logLASSO.LATLpos.final,2))); ...
    (std(logLASSO.RATL.final)/sqrt(size(logLASSO.RATL.final,2))), (std(logLASSO.RATLant.final)/sqrt(size(logLASSO.RATLant.final,2))), (std(logLASSO.RATLpos.final)/sqrt(size(logLASSO.RATLpos.final,2))); ...
    (std(SOSLASSO.LATL.final)/sqrt(size(SOSLASSO.LATL.final,2))), (std(SOSLASSO.LATLant.final)/sqrt(size(SOSLASSO.LATLant.final,2))), (std(SOSLASSO.LATLpos.final)/sqrt(size(SOSLASSO.LATLpos.final,2))); ...
    (std(SOSLASSO.RATL.final)/sqrt(size(SOSLASSO.RATL.final,2))), (std(SOSLASSO.RATLant.final)/sqrt(size(SOSLASSO.RATLant.final,2))), (std(SOSLASSO.RATLpos.final)/sqrt(size(SOSLASSO.RATLpos.final,2)))];
% bootstrapped permutation distribution for testing difference from zero
plot.bootstrappedPerm = {logLASSO.LATL.bootstrappedPerm, logLASSO.LATLant.bootstrappedPerm, logLASSO.LATLpos.bootstrappedPerm; ...
    logLASSO.RATL.bootstrappedPerm, logLASSO.RATLant.bootstrappedPerm, logLASSO.RATLpos.bootstrappedPerm; ...
    SOSLASSO.LATL.bootstrappedPerm, SOSLASSO.LATLant.bootstrappedPerm, SOSLASSO.LATLpos.bootstrappedPerm; ...
    SOSLASSO.RATL.bootstrappedPerm, SOSLASSO.RATLant.bootstrappedPerm, SOSLASSO.RATLpos.bootstrappedPerm};
% differences between methods
plot.meanDifference = [(mean(logLASSO.LATL.final) - mean(SOSLASSO.LATL.final)), (mean(logLASSO.LATLant.final) - mean(SOSLASSO.LATLant.final)), (mean(logLASSO.LATLpos.final) - mean(SOSLASSO.LATLpos.final)); ...
    (mean(logLASSO.RATL.final) - mean(SOSLASSO.RATL.final)), (mean(logLASSO.RATLant.final) - mean(SOSLASSO.RATLant.final)), (mean(logLASSO.RATLpos.final) - mean(SOSLASSO.RATLpos.final)); ...
    (mean(SOSLASSO.LATL.final) - mean(logLASSO.LATL.final)), (mean(SOSLASSO.LATLant.final) - mean(logLASSO.LATLant.final)), (mean(SOSLASSO.LATLpos.final) - mean(logLASSO.LATLpos.final)); ...
    (mean(SOSLASSO.RATL.final) - mean(logLASSO.RATL.final)), (mean(SOSLASSO.RATLant.final) - mean(logLASSO.RATLant.final)), (mean(SOSLASSO.RATLpos.final) - mean(logLASSO.RATLpos.final))];
% bootstrapped permutation distribution for testing difference from the
% other method. N.B. These bootstrapped distributions are calculated as
% (neurally-inspired penalty - vanilla LASSO). Therefore, a positive
% difference score indicates better performance by the neurally-inspired
% penalty. To get the bootstrapped permutation distribution for testing
% whether vanilla LASSO is better, multiply this by -1. 
plot.bootstrappedPermDifference = {(classificationDifference.LATL.bootstrappedPerm)*-1, (classificationDifference.LATLant.bootstrappedPerm)*-1, (classificationDifference.LATLpos.bootstrappedPerm)*-1; ...
    (classificationDifference.RATL.bootstrappedPerm)*-1, (classificationDifference.RATLant.bootstrappedPerm)*-1, (classificationDifference.RATLpos.bootstrappedPerm)*-1; ...
    classificationDifference.LATL.bootstrappedPerm, classificationDifference.LATLant.bootstrappedPerm, classificationDifference.LATLpos.bootstrappedPerm; ...
    classificationDifference.RATL.bootstrappedPerm, classificationDifference.RATLant.bootstrappedPerm, classificationDifference.RATLpos.bootstrappedPerm};
% do stats. For each row and column 
for r = 1:size(plot.mean,1)
    for c = 1:size(plot.mean,2)
        % difference from zero:
        % find the number of values in the bootstrapped permutation
        % distribution greater than the true mean
        b = sum(plot.bootstrappedPerm{r,c} > plot.mean(r,c));
        % calculate p-value using formula from Cox et al. (2024)
        plot.significance(r,c) = (b + 1)/10001;
        % difference between methods
        b = sum(plot.bootstrappedPermDifference{r,c} > plot.meanDifference(r,c));
        plot.significantDifference(r,c) = (b + 1)/10001;
    end
    % control the false discovery rate (Benjamini & Hochberg, 1995)
    plot.significance(r,:) = mafdr(plot.significance(r,:),'BHFDR',1);
    plot.significantDifference(r,:) = mafdr(plot.significantDifference(r,:),'BHFDR',1);
end

% do plots

% LASSO left
analysisIndex = 1;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar(plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_left(1:3:7,:);
% set axes
xticks([])
yticks([])
% ylabel('Accuracy')
set(gca,'FontSize',20)
ylim([0.4,1])
% plot line at chance
yline(0.5,'--','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar(1:size(plot.mean,2),plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.025;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text((1:size(plot.mean,2))-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/logLASSO_left.svg','-dsvg')
close(fig)

% LASSO right
analysisIndex = 2;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar(plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_right(1:3:7,:);
% set axes
xticks([])
yticks([])
% ylabel('Accuracy')
set(gca,'FontSize',20)
ylim([0.4,1])
% plot line at chance
yline(0.5,'--','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar(1:size(plot.mean,2),plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.025;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
t = text((1:size(plot.mean,2))-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/logLASSO_right.svg','-dsvg')
close(fig)

% SOSLASSO left
analysisIndex = 3;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar(plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_left(1:3:7,:);
% set axes
xticks([])
yticks([])
% ylabel('Accuracy')
set(gca,'FontSize',20)
ylim([0.4,1])
% plot line at chance
yline(0.5,'--','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar(1:size(plot.mean,2),plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.025;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
t = text((1:size(plot.mean,2))-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/SOSLASSO_left.svg','-dsvg')
close(fig)

% SOSLASSO right
analysisIndex = 4;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar(plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_right(1:3:7,:);
% set axes
xticks([])
yticks([])
% ylabel('Accuracy')
set(gca,'FontSize',20)
ylim([0.4,1])
% plot line at chance
yline(0.5,'--','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar(1:size(plot.mean,2),plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.025;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text((1:size(plot.mean,2))-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/SOSLASSO_right.svg','-dsvg')
close(fig)

% plot RSL. Linear LASSO and grOWL for the ATL
clear plot
% prepare the variables - mean over participants
plot.mean = [mean(linearLASSO.LATL.D1.all.final), mean(linearLASSO.LATL.D2.all.final), mean(linearLASSO.LATL.D3.all.final), mean(linearLASSO.LATLant.D1.all.final), mean(linearLASSO.LATLant.D2.all.final), mean(linearLASSO.LATLant.D3.all.final), mean(linearLASSO.LATLpos.D1.all.final), mean(linearLASSO.LATLpos.D2.all.final), mean(linearLASSO.LATLpos.D3.all.final); ...
    mean(linearLASSO.RATL.D1.all.final), mean(linearLASSO.RATL.D2.all.final), mean(linearLASSO.RATL.D3.all.final), mean(linearLASSO.RATLant.D1.all.final), mean(linearLASSO.RATLant.D2.all.final), mean(linearLASSO.RATLant.D3.all.final), mean(linearLASSO.RATLpos.D1.all.final), mean(linearLASSO.RATLpos.D2.all.final), mean(linearLASSO.RATLpos.D3.all.final); ...
    mean(grOWL.LATL.D1.all.final), mean(grOWL.LATL.D2.all.final), mean(grOWL.LATL.D3.all.final), mean(grOWL.LATLant.D1.all.final), mean(grOWL.LATLant.D2.all.final), mean(grOWL.LATLant.D3.all.final), mean(grOWL.LATLpos.D1.all.final), mean(grOWL.LATLpos.D2.all.final), mean(grOWL.LATLpos.D3.all.final); ...
    mean(grOWL.RATL.D1.all.final), mean(grOWL.RATL.D2.all.final), mean(grOWL.RATL.D3.all.final), mean(grOWL.RATLant.D1.all.final), mean(grOWL.RATLant.D2.all.final), mean(grOWL.RATLant.D3.all.final), mean(grOWL.RATLpos.D1.all.final), mean(grOWL.RATLpos.D2.all.final), mean(grOWL.RATLpos.D3.all.final)];
% standard error
plot.standardError = [(std(linearLASSO.LATL.D1.all.final)/sqrt(size(linearLASSO.LATL.D1.all.final,2))), (std(linearLASSO.LATL.D2.all.final)/sqrt(size(linearLASSO.LATL.D2.all.final,2))), (std(linearLASSO.LATL.D3.all.final)/sqrt(size(linearLASSO.LATL.D3.all.final,2))), (std(linearLASSO.LATLant.D1.all.final)/sqrt(size(linearLASSO.LATLant.D1.all.final,2))), (std(linearLASSO.LATLant.D2.all.final)/sqrt(size(linearLASSO.LATLant.D2.all.final,2))), (std(linearLASSO.LATLant.D3.all.final)/sqrt(size(linearLASSO.LATLant.D3.all.final,2))), (std(linearLASSO.LATLpos.D1.all.final)/sqrt(size(linearLASSO.LATLpos.D1.all.final,2))), (std(linearLASSO.LATLpos.D2.all.final)/sqrt(size(linearLASSO.LATLpos.D2.all.final,2))), (std(linearLASSO.LATLpos.D3.all.final)/sqrt(size(linearLASSO.LATLpos.D3.all.final,2))); 
    (std(linearLASSO.RATL.D1.all.final)/sqrt(size(linearLASSO.RATL.D1.all.final,2))), (std(linearLASSO.RATL.D2.all.final)/sqrt(size(linearLASSO.RATL.D2.all.final,2))), (std(linearLASSO.RATL.D3.all.final)/sqrt(size(linearLASSO.RATL.D3.all.final,2))), (std(linearLASSO.RATLant.D1.all.final)/sqrt(size(linearLASSO.RATLant.D1.all.final,2))), (std(linearLASSO.RATLant.D2.all.final)/sqrt(size(linearLASSO.RATLant.D2.all.final,2))), (std(linearLASSO.RATLant.D3.all.final)/sqrt(size(linearLASSO.RATLant.D3.all.final,2))), (std(linearLASSO.RATLpos.D1.all.final)/sqrt(size(linearLASSO.RATLpos.D1.all.final,2))), (std(linearLASSO.RATLpos.D2.all.final)/sqrt(size(linearLASSO.RATLpos.D2.all.final,2))), (std(linearLASSO.RATLpos.D3.all.final)/sqrt(size(linearLASSO.RATLpos.D3.all.final,2))); 
    (std(grOWL.LATL.D1.all.final)/sqrt(size(grOWL.LATL.D1.all.final,2))), (std(grOWL.LATL.D2.all.final)/sqrt(size(grOWL.LATL.D2.all.final,2))), (std(grOWL.LATL.D3.all.final)/sqrt(size(grOWL.LATL.D3.all.final,2))), (std(grOWL.LATLant.D1.all.final)/sqrt(size(grOWL.LATLant.D1.all.final,2))), (std(grOWL.LATLant.D2.all.final)/sqrt(size(grOWL.LATLant.D2.all.final,2))), (std(grOWL.LATLant.D3.all.final)/sqrt(size(grOWL.LATLant.D3.all.final,2))), (std(grOWL.LATLpos.D1.all.final)/sqrt(size(grOWL.LATLpos.D1.all.final,2))), (std(grOWL.LATLpos.D2.all.final)/sqrt(size(grOWL.LATLpos.D2.all.final,2))), (std(grOWL.LATLpos.D3.all.final)/sqrt(size(grOWL.LATLpos.D3.all.final,2))); 
    (std(grOWL.RATL.D1.all.final)/sqrt(size(grOWL.RATL.D1.all.final,2))), (std(grOWL.RATL.D2.all.final)/sqrt(size(grOWL.RATL.D2.all.final,2))), (std(grOWL.RATL.D3.all.final)/sqrt(size(grOWL.RATL.D3.all.final,2))), (std(grOWL.RATLant.D1.all.final)/sqrt(size(grOWL.RATLant.D1.all.final,2))), (std(grOWL.RATLant.D2.all.final)/sqrt(size(grOWL.RATLant.D2.all.final,2))), (std(grOWL.RATLant.D3.all.final)/sqrt(size(grOWL.RATLant.D3.all.final,2))), (std(grOWL.RATLpos.D1.all.final)/sqrt(size(grOWL.RATLpos.D1.all.final,2))), (std(grOWL.RATLpos.D2.all.final)/sqrt(size(grOWL.RATLpos.D2.all.final,2))), (std(grOWL.RATLpos.D3.all.final)/sqrt(size(grOWL.RATLpos.D3.all.final,2)))]; 
% bootstrapped permutation distribution for testing difference from zero
plot.bootstrappedPerm = {linearLASSO.LATL.D1.all.bootstrappedPerm, linearLASSO.LATL.D2.all.bootstrappedPerm,linearLASSO.LATL.D3.all.bootstrappedPerm, linearLASSO.LATLant.D1.all.bootstrappedPerm, linearLASSO.LATLant.D2.all.bootstrappedPerm,linearLASSO.LATLant.D3.all.bootstrappedPerm, linearLASSO.LATLpos.D1.all.bootstrappedPerm, linearLASSO.LATLpos.D2.all.bootstrappedPerm,linearLASSO.LATLpos.D3.all.bootstrappedPerm;
    linearLASSO.RATL.D1.all.bootstrappedPerm, linearLASSO.RATL.D2.all.bootstrappedPerm,linearLASSO.RATL.D3.all.bootstrappedPerm, linearLASSO.RATLant.D1.all.bootstrappedPerm, linearLASSO.RATLant.D2.all.bootstrappedPerm,linearLASSO.RATLant.D3.all.bootstrappedPerm, linearLASSO.RATLpos.D1.all.bootstrappedPerm, linearLASSO.RATLpos.D2.all.bootstrappedPerm,linearLASSO.RATLpos.D3.all.bootstrappedPerm;
    grOWL.LATL.D1.all.bootstrappedPerm, grOWL.LATL.D2.all.bootstrappedPerm,grOWL.LATL.D3.all.bootstrappedPerm, grOWL.LATLant.D1.all.bootstrappedPerm, grOWL.LATLant.D2.all.bootstrappedPerm,grOWL.LATLant.D3.all.bootstrappedPerm, grOWL.LATLpos.D1.all.bootstrappedPerm, grOWL.LATLpos.D2.all.bootstrappedPerm,grOWL.LATLpos.D3.all.bootstrappedPerm;
    grOWL.RATL.D1.all.bootstrappedPerm, grOWL.RATL.D2.all.bootstrappedPerm,grOWL.RATL.D3.all.bootstrappedPerm, grOWL.RATLant.D1.all.bootstrappedPerm, grOWL.RATLant.D2.all.bootstrappedPerm,grOWL.RATLant.D3.all.bootstrappedPerm, grOWL.RATLpos.D1.all.bootstrappedPerm, grOWL.RATLpos.D2.all.bootstrappedPerm,grOWL.RATLpos.D3.all.bootstrappedPerm};
% differences between methods
plot.meanDifference = [(mean(linearLASSO.LATL.D1.all.final) - mean(grOWL.LATL.D1.all.final)), (mean(linearLASSO.LATL.D2.all.final) - mean(grOWL.LATL.D2.all.final)), (mean(linearLASSO.LATL.D3.all.final) - mean(grOWL.LATL.D3.all.final)), (mean(linearLASSO.LATLant.D1.all.final) - mean(grOWL.LATLant.D1.all.final)), (mean(linearLASSO.LATLant.D2.all.final) - mean(grOWL.LATLant.D2.all.final)), (mean(linearLASSO.LATLant.D3.all.final) - mean(grOWL.LATLant.D3.all.final)), (mean(linearLASSO.LATLpos.D1.all.final) - mean(grOWL.LATLpos.D1.all.final)), (mean(linearLASSO.LATLpos.D2.all.final) - mean(grOWL.LATLpos.D2.all.final)), (mean(linearLASSO.LATLpos.D3.all.final) - mean(grOWL.LATLpos.D3.all.final)); ...
    (mean(linearLASSO.RATL.D1.all.final) - mean(grOWL.RATL.D1.all.final)), (mean(linearLASSO.RATL.D2.all.final) - mean(grOWL.RATL.D2.all.final)), (mean(linearLASSO.RATL.D3.all.final) - mean(grOWL.RATL.D3.all.final)), (mean(linearLASSO.RATLant.D1.all.final) - mean(grOWL.RATLant.D1.all.final)), (mean(linearLASSO.RATLant.D2.all.final) - mean(grOWL.RATLant.D2.all.final)), (mean(linearLASSO.RATLant.D3.all.final) - mean(grOWL.RATLant.D3.all.final)), (mean(linearLASSO.RATLpos.D1.all.final) - mean(grOWL.RATLpos.D1.all.final)), (mean(linearLASSO.RATLpos.D2.all.final) - mean(grOWL.RATLpos.D2.all.final)), (mean(linearLASSO.RATLpos.D3.all.final) - mean(grOWL.RATLpos.D3.all.final)); ...
    (mean(grOWL.LATL.D1.all.final) - mean(linearLASSO.LATL.D1.all.final)), (mean(grOWL.LATL.D2.all.final) - mean(linearLASSO.LATL.D2.all.final)), (mean(grOWL.LATL.D3.all.final) - mean(linearLASSO.LATL.D3.all.final)), (mean(grOWL.LATLant.D1.all.final) - mean(linearLASSO.LATLant.D1.all.final)), (mean(grOWL.LATLant.D2.all.final) - mean(linearLASSO.LATLant.D2.all.final)), (mean(grOWL.LATLant.D3.all.final) - mean(linearLASSO.LATLant.D3.all.final)), (mean(grOWL.LATLpos.D1.all.final) - mean(linearLASSO.LATLpos.D1.all.final)), (mean(grOWL.LATLpos.D2.all.final) - mean(linearLASSO.LATLpos.D2.all.final)), (mean(grOWL.LATLpos.D3.all.final) - mean(linearLASSO.LATLpos.D3.all.final)); ...
    (mean(grOWL.RATL.D1.all.final) - mean(linearLASSO.RATL.D1.all.final)), (mean(grOWL.RATL.D2.all.final) - mean(linearLASSO.RATL.D2.all.final)), (mean(grOWL.RATL.D3.all.final) - mean(linearLASSO.RATL.D3.all.final)), (mean(grOWL.RATLant.D1.all.final) - mean(linearLASSO.RATLant.D1.all.final)), (mean(grOWL.RATLant.D2.all.final) - mean(linearLASSO.RATLant.D2.all.final)), (mean(grOWL.RATLant.D3.all.final) - mean(linearLASSO.RATLant.D3.all.final)), (mean(grOWL.RATLpos.D1.all.final) - mean(linearLASSO.RATLpos.D1.all.final)), (mean(grOWL.RATLpos.D2.all.final) - mean(linearLASSO.RATLpos.D2.all.final)), (mean(grOWL.RATLpos.D3.all.final) - mean(linearLASSO.RATLpos.D3.all.final))];
% bootstrapped permutation distribution for testing difference from the
% other method
plot.bootstrappedPermDifference = {(RSLDifference.LATL.D1.all.bootstrappedPerm)*-1, (RSLDifference.LATL.D2.all.bootstrappedPerm)*-1, (RSLDifference.LATL.D3.all.bootstrappedPerm)*-1, (RSLDifference.LATLant.D1.all.bootstrappedPerm)*-1, (RSLDifference.LATLant.D2.all.bootstrappedPerm)*-1, (RSLDifference.LATLant.D3.all.bootstrappedPerm)*-1, (RSLDifference.LATLpos.D1.all.bootstrappedPerm)*-1, (RSLDifference.LATLpos.D2.all.bootstrappedPerm)*-1, (RSLDifference.LATLpos.D3.all.bootstrappedPerm)*-1; ...
    (RSLDifference.RATL.D1.all.bootstrappedPerm)*-1, (RSLDifference.RATL.D2.all.bootstrappedPerm)*-1, (RSLDifference.RATL.D3.all.bootstrappedPerm)*-1, (RSLDifference.RATLant.D1.all.bootstrappedPerm)*-1, (RSLDifference.RATLant.D2.all.bootstrappedPerm)*-1, (RSLDifference.RATLant.D3.all.bootstrappedPerm)*-1, (RSLDifference.RATLpos.D1.all.bootstrappedPerm)*-1, (RSLDifference.RATLpos.D2.all.bootstrappedPerm)*-1, (RSLDifference.RATLpos.D3.all.bootstrappedPerm)*-1; ...
    RSLDifference.LATL.D1.all.bootstrappedPerm, RSLDifference.LATL.D2.all.bootstrappedPerm, RSLDifference.LATL.D3.all.bootstrappedPerm, RSLDifference.LATLant.D1.all.bootstrappedPerm, RSLDifference.LATLant.D2.all.bootstrappedPerm, RSLDifference.LATLant.D3.all.bootstrappedPerm, RSLDifference.LATLpos.D1.all.bootstrappedPerm, RSLDifference.LATLpos.D2.all.bootstrappedPerm, RSLDifference.LATLpos.D3.all.bootstrappedPerm; ...
    RSLDifference.RATL.D1.all.bootstrappedPerm, RSLDifference.RATL.D2.all.bootstrappedPerm, RSLDifference.RATL.D3.all.bootstrappedPerm, RSLDifference.RATLant.D1.all.bootstrappedPerm, RSLDifference.RATLant.D2.all.bootstrappedPerm, RSLDifference.RATLant.D3.all.bootstrappedPerm, RSLDifference.RATLpos.D1.all.bootstrappedPerm, RSLDifference.RATLpos.D2.all.bootstrappedPerm, RSLDifference.RATLpos.D3.all.bootstrappedPerm};
% do stats. For each row and column 
for r = 1:size(plot.mean,1)
    for c = 1:size(plot.mean,2)
        % difference from zero:
        % find the number of values in the bootstrapped permutation
        % distribution greater than the true mean
        b = sum(plot.bootstrappedPerm{r,c} > plot.mean(r,c));
        % calculate p-value using formula from Cox et al. (2024)
        plot.significance(r,c) = (b + 1)/10001;
        % difference between methods
        b = sum(plot.bootstrappedPermDifference{r,c} > plot.meanDifference(r,c));
        plot.significantDifference(r,c) = (b + 1)/10001;
    end
    % control the false discovery rate (Benjamini & Hochberg, 1995)
    plot.significance(r,:) = mafdr(plot.significance(r,:),'BHFDR',1);
    plot.significantDifference(r,:) = mafdr(plot.significantDifference(r,:),'BHFDR',1);
end

% do plots

% LASSO left
analysisIndex = 1;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_left;
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.05;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text([1,2,3,5,6,7,9,10,11]-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/linearLASSO_left_all.svg','-dsvg')
close(fig)

% LASSO right
analysisIndex = 2;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_right;
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.05;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text([1,2,3,5,6,7,9,10,11]-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/linearLASSO_right_all.svg','-dsvg')
close(fig)

% grOWL left
analysisIndex = 3;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_left;
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.05;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text([1,2,3,5,6,7,9,10,11]-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/grOWL_left_all.svg','-dsvg')
close(fig)

% grOWL right
analysisIndex = 4;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_right;
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.05;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text([1,2,3,5,6,7,9,10,11]-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/grOWL_right_all.svg','-dsvg')
close(fig)

% animate only
clear plot
% prepare the variables - mean over participants
plot.mean = [mean(linearLASSO.LATL.D1.animate.final), mean(linearLASSO.LATL.D2.animate.final), mean(linearLASSO.LATL.D3.animate.final), mean(linearLASSO.LATLant.D1.animate.final), mean(linearLASSO.LATLant.D2.animate.final), mean(linearLASSO.LATLant.D3.animate.final), mean(linearLASSO.LATLpos.D1.animate.final), mean(linearLASSO.LATLpos.D2.animate.final), mean(linearLASSO.LATLpos.D3.animate.final); ...
    mean(linearLASSO.RATL.D1.animate.final), mean(linearLASSO.RATL.D2.animate.final), mean(linearLASSO.RATL.D3.animate.final), mean(linearLASSO.RATLant.D1.animate.final), mean(linearLASSO.RATLant.D2.animate.final), mean(linearLASSO.RATLant.D3.animate.final), mean(linearLASSO.RATLpos.D1.animate.final), mean(linearLASSO.RATLpos.D2.animate.final), mean(linearLASSO.RATLpos.D3.animate.final); ...
    mean(grOWL.LATL.D1.animate.final), mean(grOWL.LATL.D2.animate.final), mean(grOWL.LATL.D3.animate.final), mean(grOWL.LATLant.D1.animate.final), mean(grOWL.LATLant.D2.animate.final), mean(grOWL.LATLant.D3.animate.final), mean(grOWL.LATLpos.D1.animate.final), mean(grOWL.LATLpos.D2.animate.final), mean(grOWL.LATLpos.D3.animate.final); ...
    mean(grOWL.RATL.D1.animate.final), mean(grOWL.RATL.D2.animate.final), mean(grOWL.RATL.D3.animate.final), mean(grOWL.RATLant.D1.animate.final), mean(grOWL.RATLant.D2.animate.final), mean(grOWL.RATLant.D3.animate.final), mean(grOWL.RATLpos.D1.animate.final), mean(grOWL.RATLpos.D2.animate.final), mean(grOWL.RATLpos.D3.animate.final)];
% standard error
plot.standardError = [(std(linearLASSO.LATL.D1.animate.final)/sqrt(size(linearLASSO.LATL.D1.animate.final,2))), (std(linearLASSO.LATL.D2.animate.final)/sqrt(size(linearLASSO.LATL.D2.animate.final,2))), (std(linearLASSO.LATL.D3.animate.final)/sqrt(size(linearLASSO.LATL.D3.animate.final,2))), (std(linearLASSO.LATLant.D1.animate.final)/sqrt(size(linearLASSO.LATLant.D1.animate.final,2))), (std(linearLASSO.LATLant.D2.animate.final)/sqrt(size(linearLASSO.LATLant.D2.animate.final,2))), (std(linearLASSO.LATLant.D3.animate.final)/sqrt(size(linearLASSO.LATLant.D3.animate.final,2))), (std(linearLASSO.LATLpos.D1.animate.final)/sqrt(size(linearLASSO.LATLpos.D1.animate.final,2))), (std(linearLASSO.LATLpos.D2.animate.final)/sqrt(size(linearLASSO.LATLpos.D2.animate.final,2))), (std(linearLASSO.LATLpos.D3.animate.final)/sqrt(size(linearLASSO.LATLpos.D3.animate.final,2))); 
    (std(linearLASSO.RATL.D1.animate.final)/sqrt(size(linearLASSO.RATL.D1.animate.final,2))), (std(linearLASSO.RATL.D2.animate.final)/sqrt(size(linearLASSO.RATL.D2.animate.final,2))), (std(linearLASSO.RATL.D3.animate.final)/sqrt(size(linearLASSO.RATL.D3.animate.final,2))), (std(linearLASSO.RATLant.D1.animate.final)/sqrt(size(linearLASSO.RATLant.D1.animate.final,2))), (std(linearLASSO.RATLant.D2.animate.final)/sqrt(size(linearLASSO.RATLant.D2.animate.final,2))), (std(linearLASSO.RATLant.D3.animate.final)/sqrt(size(linearLASSO.RATLant.D3.animate.final,2))), (std(linearLASSO.RATLpos.D1.animate.final)/sqrt(size(linearLASSO.RATLpos.D1.animate.final,2))), (std(linearLASSO.RATLpos.D2.animate.final)/sqrt(size(linearLASSO.RATLpos.D2.animate.final,2))), (std(linearLASSO.RATLpos.D3.animate.final)/sqrt(size(linearLASSO.RATLpos.D3.animate.final,2))); 
    (std(grOWL.LATL.D1.animate.final)/sqrt(size(grOWL.LATL.D1.animate.final,2))), (std(grOWL.LATL.D2.animate.final)/sqrt(size(grOWL.LATL.D2.animate.final,2))), (std(grOWL.LATL.D3.animate.final)/sqrt(size(grOWL.LATL.D3.animate.final,2))), (std(grOWL.LATLant.D1.animate.final)/sqrt(size(grOWL.LATLant.D1.animate.final,2))), (std(grOWL.LATLant.D2.animate.final)/sqrt(size(grOWL.LATLant.D2.animate.final,2))), (std(grOWL.LATLant.D3.animate.final)/sqrt(size(grOWL.LATLant.D3.animate.final,2))), (std(grOWL.LATLpos.D1.animate.final)/sqrt(size(grOWL.LATLpos.D1.animate.final,2))), (std(grOWL.LATLpos.D2.animate.final)/sqrt(size(grOWL.LATLpos.D2.animate.final,2))), (std(grOWL.LATLpos.D3.animate.final)/sqrt(size(grOWL.LATLpos.D3.animate.final,2))); 
    (std(grOWL.RATL.D1.animate.final)/sqrt(size(grOWL.RATL.D1.animate.final,2))), (std(grOWL.RATL.D2.animate.final)/sqrt(size(grOWL.RATL.D2.animate.final,2))), (std(grOWL.RATL.D3.animate.final)/sqrt(size(grOWL.RATL.D3.animate.final,2))), (std(grOWL.RATLant.D1.animate.final)/sqrt(size(grOWL.RATLant.D1.animate.final,2))), (std(grOWL.RATLant.D2.animate.final)/sqrt(size(grOWL.RATLant.D2.animate.final,2))), (std(grOWL.RATLant.D3.animate.final)/sqrt(size(grOWL.RATLant.D3.animate.final,2))), (std(grOWL.RATLpos.D1.animate.final)/sqrt(size(grOWL.RATLpos.D1.animate.final,2))), (std(grOWL.RATLpos.D2.animate.final)/sqrt(size(grOWL.RATLpos.D2.animate.final,2))), (std(grOWL.RATLpos.D3.animate.final)/sqrt(size(grOWL.RATLpos.D3.animate.final,2)))]; 
% bootstrapped permutation distribution for testing difference from zero
plot.bootstrappedPerm = {linearLASSO.LATL.D1.animate.bootstrappedPerm, linearLASSO.LATL.D2.animate.bootstrappedPerm,linearLASSO.LATL.D3.animate.bootstrappedPerm, linearLASSO.LATLant.D1.animate.bootstrappedPerm, linearLASSO.LATLant.D2.animate.bootstrappedPerm,linearLASSO.LATLant.D3.animate.bootstrappedPerm, linearLASSO.LATLpos.D1.animate.bootstrappedPerm, linearLASSO.LATLpos.D2.animate.bootstrappedPerm,linearLASSO.LATLpos.D3.animate.bootstrappedPerm;
    linearLASSO.RATL.D1.animate.bootstrappedPerm, linearLASSO.RATL.D2.animate.bootstrappedPerm,linearLASSO.RATL.D3.animate.bootstrappedPerm, linearLASSO.RATLant.D1.animate.bootstrappedPerm, linearLASSO.RATLant.D2.animate.bootstrappedPerm,linearLASSO.RATLant.D3.animate.bootstrappedPerm, linearLASSO.RATLpos.D1.animate.bootstrappedPerm, linearLASSO.RATLpos.D2.animate.bootstrappedPerm,linearLASSO.RATLpos.D3.animate.bootstrappedPerm;
    grOWL.LATL.D1.animate.bootstrappedPerm, grOWL.LATL.D2.animate.bootstrappedPerm,grOWL.LATL.D3.animate.bootstrappedPerm, grOWL.LATLant.D1.animate.bootstrappedPerm, grOWL.LATLant.D2.animate.bootstrappedPerm,grOWL.LATLant.D3.animate.bootstrappedPerm, grOWL.LATLpos.D1.animate.bootstrappedPerm, grOWL.LATLpos.D2.animate.bootstrappedPerm,grOWL.LATLpos.D3.animate.bootstrappedPerm;
    grOWL.RATL.D1.animate.bootstrappedPerm, grOWL.RATL.D2.animate.bootstrappedPerm,grOWL.RATL.D3.animate.bootstrappedPerm, grOWL.RATLant.D1.animate.bootstrappedPerm, grOWL.RATLant.D2.animate.bootstrappedPerm,grOWL.RATLant.D3.animate.bootstrappedPerm, grOWL.RATLpos.D1.animate.bootstrappedPerm, grOWL.RATLpos.D2.animate.bootstrappedPerm,grOWL.RATLpos.D3.animate.bootstrappedPerm};
% differences between methods
plot.meanDifference = [(mean(linearLASSO.LATL.D1.animate.final) - mean(grOWL.LATL.D1.animate.final)), (mean(linearLASSO.LATL.D2.animate.final) - mean(grOWL.LATL.D2.animate.final)), (mean(linearLASSO.LATL.D3.animate.final) - mean(grOWL.LATL.D3.animate.final)), (mean(linearLASSO.LATLant.D1.animate.final) - mean(grOWL.LATLant.D1.animate.final)), (mean(linearLASSO.LATLant.D2.animate.final) - mean(grOWL.LATLant.D2.animate.final)), (mean(linearLASSO.LATLant.D3.animate.final) - mean(grOWL.LATLant.D3.animate.final)), (mean(linearLASSO.LATLpos.D1.animate.final) - mean(grOWL.LATLpos.D1.animate.final)), (mean(linearLASSO.LATLpos.D2.animate.final) - mean(grOWL.LATLpos.D2.animate.final)), (mean(linearLASSO.LATLpos.D3.animate.final) - mean(grOWL.LATLpos.D3.animate.final)); ...
    (mean(linearLASSO.RATL.D1.animate.final) - mean(grOWL.RATL.D1.animate.final)), (mean(linearLASSO.RATL.D2.animate.final) - mean(grOWL.RATL.D2.animate.final)), (mean(linearLASSO.RATL.D3.animate.final) - mean(grOWL.RATL.D3.animate.final)), (mean(linearLASSO.RATLant.D1.animate.final) - mean(grOWL.RATLant.D1.animate.final)), (mean(linearLASSO.RATLant.D2.animate.final) - mean(grOWL.RATLant.D2.animate.final)), (mean(linearLASSO.RATLant.D3.animate.final) - mean(grOWL.RATLant.D3.animate.final)), (mean(linearLASSO.RATLpos.D1.animate.final) - mean(grOWL.RATLpos.D1.animate.final)), (mean(linearLASSO.RATLpos.D2.animate.final) - mean(grOWL.RATLpos.D2.animate.final)), (mean(linearLASSO.RATLpos.D3.animate.final) - mean(grOWL.RATLpos.D3.animate.final)); ...
    (mean(grOWL.LATL.D1.animate.final) - mean(linearLASSO.LATL.D1.animate.final)), (mean(grOWL.LATL.D2.animate.final) - mean(linearLASSO.LATL.D2.animate.final)), (mean(grOWL.LATL.D3.animate.final) - mean(linearLASSO.LATL.D3.animate.final)), (mean(grOWL.LATLant.D1.animate.final) - mean(linearLASSO.LATLant.D1.animate.final)), (mean(grOWL.LATLant.D2.animate.final) - mean(linearLASSO.LATLant.D2.animate.final)), (mean(grOWL.LATLant.D3.animate.final) - mean(linearLASSO.LATLant.D3.animate.final)), (mean(grOWL.LATLpos.D1.animate.final) - mean(linearLASSO.LATLpos.D1.animate.final)), (mean(grOWL.LATLpos.D2.animate.final) - mean(linearLASSO.LATLpos.D2.animate.final)), (mean(grOWL.LATLpos.D3.animate.final) - mean(linearLASSO.LATLpos.D3.animate.final)); ...
    (mean(grOWL.RATL.D1.animate.final) - mean(linearLASSO.RATL.D1.animate.final)), (mean(grOWL.RATL.D2.animate.final) - mean(linearLASSO.RATL.D2.animate.final)), (mean(grOWL.RATL.D3.animate.final) - mean(linearLASSO.RATL.D3.animate.final)), (mean(grOWL.RATLant.D1.animate.final) - mean(linearLASSO.RATLant.D1.animate.final)), (mean(grOWL.RATLant.D2.animate.final) - mean(linearLASSO.RATLant.D2.animate.final)), (mean(grOWL.RATLant.D3.animate.final) - mean(linearLASSO.RATLant.D3.animate.final)), (mean(grOWL.RATLpos.D1.animate.final) - mean(linearLASSO.RATLpos.D1.animate.final)), (mean(grOWL.RATLpos.D2.animate.final) - mean(linearLASSO.RATLpos.D2.animate.final)), (mean(grOWL.RATLpos.D3.animate.final) - mean(linearLASSO.RATLpos.D3.animate.final))];
% bootstrapped permutation distribution for testing difference from the
% other method
plot.bootstrappedPermDifference = {(RSLDifference.LATL.D1.animate.bootstrappedPerm)*-1, (RSLDifference.LATL.D2.animate.bootstrappedPerm)*-1, (RSLDifference.LATL.D3.animate.bootstrappedPerm)*-1, (RSLDifference.LATLant.D1.animate.bootstrappedPerm)*-1, (RSLDifference.LATLant.D2.animate.bootstrappedPerm)*-1, (RSLDifference.LATLant.D3.animate.bootstrappedPerm)*-1, (RSLDifference.LATLpos.D1.animate.bootstrappedPerm)*-1, (RSLDifference.LATLpos.D2.animate.bootstrappedPerm)*-1, (RSLDifference.LATLpos.D3.animate.bootstrappedPerm)*-1; ...
    (RSLDifference.RATL.D1.animate.bootstrappedPerm)*-1, (RSLDifference.RATL.D2.animate.bootstrappedPerm)*-1, (RSLDifference.RATL.D3.animate.bootstrappedPerm)*-1, (RSLDifference.RATLant.D1.animate.bootstrappedPerm)*-1, (RSLDifference.RATLant.D2.animate.bootstrappedPerm)*-1, (RSLDifference.RATLant.D3.animate.bootstrappedPerm)*-1, (RSLDifference.RATLpos.D1.animate.bootstrappedPerm)*-1, (RSLDifference.RATLpos.D2.animate.bootstrappedPerm)*-1, (RSLDifference.RATLpos.D3.animate.bootstrappedPerm)*-1; ...
    RSLDifference.LATL.D1.animate.bootstrappedPerm, RSLDifference.LATL.D2.animate.bootstrappedPerm, RSLDifference.LATL.D3.animate.bootstrappedPerm, RSLDifference.LATLant.D1.animate.bootstrappedPerm, RSLDifference.LATLant.D2.animate.bootstrappedPerm, RSLDifference.LATLant.D3.animate.bootstrappedPerm, RSLDifference.LATLpos.D1.animate.bootstrappedPerm, RSLDifference.LATLpos.D2.animate.bootstrappedPerm, RSLDifference.LATLpos.D3.animate.bootstrappedPerm; ...
    RSLDifference.RATL.D1.animate.bootstrappedPerm, RSLDifference.RATL.D2.animate.bootstrappedPerm, RSLDifference.RATL.D3.animate.bootstrappedPerm, RSLDifference.RATLant.D1.animate.bootstrappedPerm, RSLDifference.RATLant.D2.animate.bootstrappedPerm, RSLDifference.RATLant.D3.animate.bootstrappedPerm, RSLDifference.RATLpos.D1.animate.bootstrappedPerm, RSLDifference.RATLpos.D2.animate.bootstrappedPerm, RSLDifference.RATLpos.D3.animate.bootstrappedPerm};
% do stats. For each row and column 
for r = 1:size(plot.mean,1)
    for c = 1:size(plot.mean,2)
        % difference from zero:
        % find the number of values in the bootstrapped permutation
        % distribution greater than the true mean
        b = sum(plot.bootstrappedPerm{r,c} > plot.mean(r,c));
        % calculate p-value using formula from Cox et al. (2024)
        plot.significance(r,c) = (b + 1)/10001;
        % difference between methods
        b = sum(plot.bootstrappedPermDifference{r,c} > plot.meanDifference(r,c));
        plot.significantDifference(r,c) = (b + 1)/10001;
    end
    % control the false discovery rate (Benjamini & Hochberg, 1995)
    plot.significance(r,:) = mafdr(plot.significance(r,:),'BHFDR',1);
    plot.significantDifference(r,:) = mafdr(plot.significantDifference(r,:),'BHFDR',1);
end

% do plots

% LASSO left
analysisIndex = 1;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_left;
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,0.6])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.025;
positionOnPage(positionOnPage > 0.575) = 0.575;
% plot the stars
text([1,2,3,5,6,7,9,10,11]-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/linearLASSO_left_animate.svg','-dsvg')
close(fig)

% LASSO right
analysisIndex = 2;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_right;
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,0.6])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.025;
positionOnPage(positionOnPage > 0.575) = 0.575;
% plot the stars
text([1,2,3,5,6,7,9,10,11]-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/linearLASSO_right_animate.svg','-dsvg')
close(fig)

% grOWL left
analysisIndex = 3;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_left;
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,0.6])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.025;
positionOnPage(positionOnPage > 0.575) = 0.575;
% plot the stars
text([1,2,3,5,6,7,9,10,11]-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/grOWL_left_animate.svg','-dsvg')
close(fig)

% grOWL right
analysisIndex = 4;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_right;
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,0.6])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.025;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text([1,2,3,5,6,7,9,10,11]-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/grOWL_right_animate.svg','-dsvg')
close(fig)

% inanimate only
clear plot
% prepare the variables - mean over participants
plot.mean = [mean(linearLASSO.LATL.D1.inanimate.final), mean(linearLASSO.LATL.D2.inanimate.final), mean(linearLASSO.LATL.D3.inanimate.final), mean(linearLASSO.LATLant.D1.inanimate.final), mean(linearLASSO.LATLant.D2.inanimate.final), mean(linearLASSO.LATLant.D3.inanimate.final), mean(linearLASSO.LATLpos.D1.inanimate.final), mean(linearLASSO.LATLpos.D2.inanimate.final), mean(linearLASSO.LATLpos.D3.inanimate.final); ...
    mean(linearLASSO.RATL.D1.inanimate.final), mean(linearLASSO.RATL.D2.inanimate.final), mean(linearLASSO.RATL.D3.inanimate.final), mean(linearLASSO.RATLant.D1.inanimate.final), mean(linearLASSO.RATLant.D2.inanimate.final), mean(linearLASSO.RATLant.D3.inanimate.final), mean(linearLASSO.RATLpos.D1.inanimate.final), mean(linearLASSO.RATLpos.D2.inanimate.final), mean(linearLASSO.RATLpos.D3.inanimate.final); ...
    mean(grOWL.LATL.D1.inanimate.final), mean(grOWL.LATL.D2.inanimate.final), mean(grOWL.LATL.D3.inanimate.final), mean(grOWL.LATLant.D1.inanimate.final), mean(grOWL.LATLant.D2.inanimate.final), mean(grOWL.LATLant.D3.inanimate.final), mean(grOWL.LATLpos.D1.inanimate.final), mean(grOWL.LATLpos.D2.inanimate.final), mean(grOWL.LATLpos.D3.inanimate.final); ...
    mean(grOWL.RATL.D1.inanimate.final), mean(grOWL.RATL.D2.inanimate.final), mean(grOWL.RATL.D3.inanimate.final), mean(grOWL.RATLant.D1.inanimate.final), mean(grOWL.RATLant.D2.inanimate.final), mean(grOWL.RATLant.D3.inanimate.final), mean(grOWL.RATLpos.D1.inanimate.final), mean(grOWL.RATLpos.D2.inanimate.final), mean(grOWL.RATLpos.D3.inanimate.final)];
% standard error
plot.standardError = [(std(linearLASSO.LATL.D1.inanimate.final)/sqrt(size(linearLASSO.LATL.D1.inanimate.final,2))), (std(linearLASSO.LATL.D2.inanimate.final)/sqrt(size(linearLASSO.LATL.D2.inanimate.final,2))), (std(linearLASSO.LATL.D3.inanimate.final)/sqrt(size(linearLASSO.LATL.D3.inanimate.final,2))), (std(linearLASSO.LATLant.D1.inanimate.final)/sqrt(size(linearLASSO.LATLant.D1.inanimate.final,2))), (std(linearLASSO.LATLant.D2.inanimate.final)/sqrt(size(linearLASSO.LATLant.D2.inanimate.final,2))), (std(linearLASSO.LATLant.D3.inanimate.final)/sqrt(size(linearLASSO.LATLant.D3.inanimate.final,2))), (std(linearLASSO.LATLpos.D1.inanimate.final)/sqrt(size(linearLASSO.LATLpos.D1.inanimate.final,2))), (std(linearLASSO.LATLpos.D2.inanimate.final)/sqrt(size(linearLASSO.LATLpos.D2.inanimate.final,2))), (std(linearLASSO.LATLpos.D3.inanimate.final)/sqrt(size(linearLASSO.LATLpos.D3.inanimate.final,2))); 
    (std(linearLASSO.RATL.D1.inanimate.final)/sqrt(size(linearLASSO.RATL.D1.inanimate.final,2))), (std(linearLASSO.RATL.D2.inanimate.final)/sqrt(size(linearLASSO.RATL.D2.inanimate.final,2))), (std(linearLASSO.RATL.D3.inanimate.final)/sqrt(size(linearLASSO.RATL.D3.inanimate.final,2))), (std(linearLASSO.RATLant.D1.inanimate.final)/sqrt(size(linearLASSO.RATLant.D1.inanimate.final,2))), (std(linearLASSO.RATLant.D2.inanimate.final)/sqrt(size(linearLASSO.RATLant.D2.inanimate.final,2))), (std(linearLASSO.RATLant.D3.inanimate.final)/sqrt(size(linearLASSO.RATLant.D3.inanimate.final,2))), (std(linearLASSO.RATLpos.D1.inanimate.final)/sqrt(size(linearLASSO.RATLpos.D1.inanimate.final,2))), (std(linearLASSO.RATLpos.D2.inanimate.final)/sqrt(size(linearLASSO.RATLpos.D2.inanimate.final,2))), (std(linearLASSO.RATLpos.D3.inanimate.final)/sqrt(size(linearLASSO.RATLpos.D3.inanimate.final,2))); 
    (std(grOWL.LATL.D1.inanimate.final)/sqrt(size(grOWL.LATL.D1.inanimate.final,2))), (std(grOWL.LATL.D2.inanimate.final)/sqrt(size(grOWL.LATL.D2.inanimate.final,2))), (std(grOWL.LATL.D3.inanimate.final)/sqrt(size(grOWL.LATL.D3.inanimate.final,2))), (std(grOWL.LATLant.D1.inanimate.final)/sqrt(size(grOWL.LATLant.D1.inanimate.final,2))), (std(grOWL.LATLant.D2.inanimate.final)/sqrt(size(grOWL.LATLant.D2.inanimate.final,2))), (std(grOWL.LATLant.D3.inanimate.final)/sqrt(size(grOWL.LATLant.D3.inanimate.final,2))), (std(grOWL.LATLpos.D1.inanimate.final)/sqrt(size(grOWL.LATLpos.D1.inanimate.final,2))), (std(grOWL.LATLpos.D2.inanimate.final)/sqrt(size(grOWL.LATLpos.D2.inanimate.final,2))), (std(grOWL.LATLpos.D3.inanimate.final)/sqrt(size(grOWL.LATLpos.D3.inanimate.final,2))); 
    (std(grOWL.RATL.D1.inanimate.final)/sqrt(size(grOWL.RATL.D1.inanimate.final,2))), (std(grOWL.RATL.D2.inanimate.final)/sqrt(size(grOWL.RATL.D2.inanimate.final,2))), (std(grOWL.RATL.D3.inanimate.final)/sqrt(size(grOWL.RATL.D3.inanimate.final,2))), (std(grOWL.RATLant.D1.inanimate.final)/sqrt(size(grOWL.RATLant.D1.inanimate.final,2))), (std(grOWL.RATLant.D2.inanimate.final)/sqrt(size(grOWL.RATLant.D2.inanimate.final,2))), (std(grOWL.RATLant.D3.inanimate.final)/sqrt(size(grOWL.RATLant.D3.inanimate.final,2))), (std(grOWL.RATLpos.D1.inanimate.final)/sqrt(size(grOWL.RATLpos.D1.inanimate.final,2))), (std(grOWL.RATLpos.D2.inanimate.final)/sqrt(size(grOWL.RATLpos.D2.inanimate.final,2))), (std(grOWL.RATLpos.D3.inanimate.final)/sqrt(size(grOWL.RATLpos.D3.inanimate.final,2)))]; 
% bootstrapped permutation distribution for testing difference from zero
plot.bootstrappedPerm = {linearLASSO.LATL.D1.inanimate.bootstrappedPerm, linearLASSO.LATL.D2.inanimate.bootstrappedPerm,linearLASSO.LATL.D3.inanimate.bootstrappedPerm, linearLASSO.LATLant.D1.inanimate.bootstrappedPerm, linearLASSO.LATLant.D2.inanimate.bootstrappedPerm,linearLASSO.LATLant.D3.inanimate.bootstrappedPerm, linearLASSO.LATLpos.D1.inanimate.bootstrappedPerm, linearLASSO.LATLpos.D2.inanimate.bootstrappedPerm,linearLASSO.LATLpos.D3.inanimate.bootstrappedPerm;
    linearLASSO.RATL.D1.inanimate.bootstrappedPerm, linearLASSO.RATL.D2.inanimate.bootstrappedPerm,linearLASSO.RATL.D3.inanimate.bootstrappedPerm, linearLASSO.RATLant.D1.inanimate.bootstrappedPerm, linearLASSO.RATLant.D2.inanimate.bootstrappedPerm,linearLASSO.RATLant.D3.inanimate.bootstrappedPerm, linearLASSO.RATLpos.D1.inanimate.bootstrappedPerm, linearLASSO.RATLpos.D2.inanimate.bootstrappedPerm,linearLASSO.RATLpos.D3.inanimate.bootstrappedPerm;
    grOWL.LATL.D1.inanimate.bootstrappedPerm, grOWL.LATL.D2.inanimate.bootstrappedPerm,grOWL.LATL.D3.inanimate.bootstrappedPerm, grOWL.LATLant.D1.inanimate.bootstrappedPerm, grOWL.LATLant.D2.inanimate.bootstrappedPerm,grOWL.LATLant.D3.inanimate.bootstrappedPerm, grOWL.LATLpos.D1.inanimate.bootstrappedPerm, grOWL.LATLpos.D2.inanimate.bootstrappedPerm,grOWL.LATLpos.D3.inanimate.bootstrappedPerm;
    grOWL.RATL.D1.inanimate.bootstrappedPerm, grOWL.RATL.D2.inanimate.bootstrappedPerm,grOWL.RATL.D3.inanimate.bootstrappedPerm, grOWL.RATLant.D1.inanimate.bootstrappedPerm, grOWL.RATLant.D2.inanimate.bootstrappedPerm,grOWL.RATLant.D3.inanimate.bootstrappedPerm, grOWL.RATLpos.D1.inanimate.bootstrappedPerm, grOWL.RATLpos.D2.inanimate.bootstrappedPerm,grOWL.RATLpos.D3.inanimate.bootstrappedPerm};
% differences between methods
plot.meanDifference = [(mean(linearLASSO.LATL.D1.inanimate.final) - mean(grOWL.LATL.D1.inanimate.final)), (mean(linearLASSO.LATL.D2.inanimate.final) - mean(grOWL.LATL.D2.inanimate.final)), (mean(linearLASSO.LATL.D3.inanimate.final) - mean(grOWL.LATL.D3.inanimate.final)), (mean(linearLASSO.LATLant.D1.inanimate.final) - mean(grOWL.LATLant.D1.inanimate.final)), (mean(linearLASSO.LATLant.D2.inanimate.final) - mean(grOWL.LATLant.D2.inanimate.final)), (mean(linearLASSO.LATLant.D3.inanimate.final) - mean(grOWL.LATLant.D3.inanimate.final)), (mean(linearLASSO.LATLpos.D1.inanimate.final) - mean(grOWL.LATLpos.D1.inanimate.final)), (mean(linearLASSO.LATLpos.D2.inanimate.final) - mean(grOWL.LATLpos.D2.inanimate.final)), (mean(linearLASSO.LATLpos.D3.inanimate.final) - mean(grOWL.LATLpos.D3.inanimate.final)); ...
    (mean(linearLASSO.RATL.D1.inanimate.final) - mean(grOWL.RATL.D1.inanimate.final)), (mean(linearLASSO.RATL.D2.inanimate.final) - mean(grOWL.RATL.D2.inanimate.final)), (mean(linearLASSO.RATL.D3.inanimate.final) - mean(grOWL.RATL.D3.inanimate.final)), (mean(linearLASSO.RATLant.D1.inanimate.final) - mean(grOWL.RATLant.D1.inanimate.final)), (mean(linearLASSO.RATLant.D2.inanimate.final) - mean(grOWL.RATLant.D2.inanimate.final)), (mean(linearLASSO.RATLant.D3.inanimate.final) - mean(grOWL.RATLant.D3.inanimate.final)), (mean(linearLASSO.RATLpos.D1.inanimate.final) - mean(grOWL.RATLpos.D1.inanimate.final)), (mean(linearLASSO.RATLpos.D2.inanimate.final) - mean(grOWL.RATLpos.D2.inanimate.final)), (mean(linearLASSO.RATLpos.D3.inanimate.final) - mean(grOWL.RATLpos.D3.inanimate.final)); ...
    (mean(grOWL.LATL.D1.inanimate.final) - mean(linearLASSO.LATL.D1.inanimate.final)), (mean(grOWL.LATL.D2.inanimate.final) - mean(linearLASSO.LATL.D2.inanimate.final)), (mean(grOWL.LATL.D3.inanimate.final) - mean(linearLASSO.LATL.D3.inanimate.final)), (mean(grOWL.LATLant.D1.inanimate.final) - mean(linearLASSO.LATLant.D1.inanimate.final)), (mean(grOWL.LATLant.D2.inanimate.final) - mean(linearLASSO.LATLant.D2.inanimate.final)), (mean(grOWL.LATLant.D3.inanimate.final) - mean(linearLASSO.LATLant.D3.inanimate.final)), (mean(grOWL.LATLpos.D1.inanimate.final) - mean(linearLASSO.LATLpos.D1.inanimate.final)), (mean(grOWL.LATLpos.D2.inanimate.final) - mean(linearLASSO.LATLpos.D2.inanimate.final)), (mean(grOWL.LATLpos.D3.inanimate.final) - mean(linearLASSO.LATLpos.D3.inanimate.final)); ...
    (mean(grOWL.RATL.D1.inanimate.final) - mean(linearLASSO.RATL.D1.inanimate.final)), (mean(grOWL.RATL.D2.inanimate.final) - mean(linearLASSO.RATL.D2.inanimate.final)), (mean(grOWL.RATL.D3.inanimate.final) - mean(linearLASSO.RATL.D3.inanimate.final)), (mean(grOWL.RATLant.D1.inanimate.final) - mean(linearLASSO.RATLant.D1.inanimate.final)), (mean(grOWL.RATLant.D2.inanimate.final) - mean(linearLASSO.RATLant.D2.inanimate.final)), (mean(grOWL.RATLant.D3.inanimate.final) - mean(linearLASSO.RATLant.D3.inanimate.final)), (mean(grOWL.RATLpos.D1.inanimate.final) - mean(linearLASSO.RATLpos.D1.inanimate.final)), (mean(grOWL.RATLpos.D2.inanimate.final) - mean(linearLASSO.RATLpos.D2.inanimate.final)), (mean(grOWL.RATLpos.D3.inanimate.final) - mean(linearLASSO.RATLpos.D3.inanimate.final))];
% bootstrapped permutation distribution for testing difference from the
% other method
plot.bootstrappedPermDifference = {(RSLDifference.LATL.D1.inanimate.bootstrappedPerm)*-1, (RSLDifference.LATL.D2.inanimate.bootstrappedPerm)*-1, (RSLDifference.LATL.D3.inanimate.bootstrappedPerm)*-1, (RSLDifference.LATLant.D1.inanimate.bootstrappedPerm)*-1, (RSLDifference.LATLant.D2.inanimate.bootstrappedPerm)*-1, (RSLDifference.LATLant.D3.inanimate.bootstrappedPerm)*-1, (RSLDifference.LATLpos.D1.inanimate.bootstrappedPerm)*-1, (RSLDifference.LATLpos.D2.inanimate.bootstrappedPerm)*-1, (RSLDifference.LATLpos.D3.inanimate.bootstrappedPerm)*-1; ...
    (RSLDifference.RATL.D1.inanimate.bootstrappedPerm)*-1, (RSLDifference.RATL.D2.inanimate.bootstrappedPerm)*-1, (RSLDifference.RATL.D3.inanimate.bootstrappedPerm)*-1, (RSLDifference.RATLant.D1.inanimate.bootstrappedPerm)*-1, (RSLDifference.RATLant.D2.inanimate.bootstrappedPerm)*-1, (RSLDifference.RATLant.D3.inanimate.bootstrappedPerm)*-1, (RSLDifference.RATLpos.D1.inanimate.bootstrappedPerm)*-1, (RSLDifference.RATLpos.D2.inanimate.bootstrappedPerm)*-1, (RSLDifference.RATLpos.D3.inanimate.bootstrappedPerm)*-1; ...
    RSLDifference.LATL.D1.inanimate.bootstrappedPerm, RSLDifference.LATL.D2.inanimate.bootstrappedPerm, RSLDifference.LATL.D3.inanimate.bootstrappedPerm, RSLDifference.LATLant.D1.inanimate.bootstrappedPerm, RSLDifference.LATLant.D2.inanimate.bootstrappedPerm, RSLDifference.LATLant.D3.inanimate.bootstrappedPerm, RSLDifference.LATLpos.D1.inanimate.bootstrappedPerm, RSLDifference.LATLpos.D2.inanimate.bootstrappedPerm, RSLDifference.LATLpos.D3.inanimate.bootstrappedPerm; ...
    RSLDifference.RATL.D1.inanimate.bootstrappedPerm, RSLDifference.RATL.D2.inanimate.bootstrappedPerm, RSLDifference.RATL.D3.inanimate.bootstrappedPerm, RSLDifference.RATLant.D1.inanimate.bootstrappedPerm, RSLDifference.RATLant.D2.inanimate.bootstrappedPerm, RSLDifference.RATLant.D3.inanimate.bootstrappedPerm, RSLDifference.RATLpos.D1.inanimate.bootstrappedPerm, RSLDifference.RATLpos.D2.inanimate.bootstrappedPerm, RSLDifference.RATLpos.D3.inanimate.bootstrappedPerm};
% do stats. For each row and column 
for r = 1:size(plot.mean,1)
    for c = 1:size(plot.mean,2)
        % difference from zero:
        % find the number of values in the bootstrapped permutation
        % distribution greater than the true mean
        b = sum(plot.bootstrappedPerm{r,c} > plot.mean(r,c));
        % calculate p-value using formula from Cox et al. (2024)
        plot.significance(r,c) = (b + 1)/10001;
        % difference between methods
        b = sum(plot.bootstrappedPermDifference{r,c} > plot.meanDifference(r,c));
        plot.significantDifference(r,c) = (b + 1)/10001;
    end
    % control the false discovery rate (Benjamini & Hochberg, 1995)
    plot.significance(r,:) = mafdr(plot.significance(r,:),'BHFDR',1);
    plot.significantDifference(r,:) = mafdr(plot.significantDifference(r,:),'BHFDR',1);
end

% do plots
% LASSO left
analysisIndex = 1;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_left;
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,0.6])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.025;
positionOnPage(positionOnPage > 0.575) = 0.575;
% plot the stars
text([1,2,3,5,6,7,9,10,11]-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/linearLASSO_left_inanimate.svg','-dsvg')
close(fig)

% LASSO right
analysisIndex = 2;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_right;
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,0.6])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.025;
positionOnPage(positionOnPage > 0.575) = 0.575;
% plot the stars
text([1,2,3,5,6,7,9,10,11]-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/linearLASSO_right_inanimate.svg','-dsvg')
close(fig)

% grOWL left
analysisIndex = 3;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_left;
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,0.6])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.025;
positionOnPage(positionOnPage > 0.575) = 0.575;
% plot the stars
text([1,2,3,5,6,7,9,10,11]-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/grOWL_left_inanimate.svg','-dsvg')
close(fig)

% grOWL right
analysisIndex = 4;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_right;
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,0.6])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar([1,2,3,5,6,7,9,10,11],plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.025;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text([1,2,3,5,6,7,9,10,11]-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/grOWL_right_inanimate.svg','-dsvg')
close(fig)

% whole-brain classification
% plot classifications. LASSO and SOSLASSO for the ATL
clear plot
% prepare the variables - mean over participants
plot.mean = [mean(logLASSO.wholebrain.final); ...
    mean(SOSLASSO.wholebrain.final)];
% standard error
plot.standardError = [(std(logLASSO.wholebrain.final)/sqrt(size(logLASSO.wholebrain.final,2))); 
    (std(SOSLASSO.wholebrain.final)/sqrt(size(SOSLASSO.RATL.final,2)))];
 % bootstrapped permutation distribution for testing difference from zero
plot.bootstrappedPerm = {logLASSO.wholebrain.bootstrappedPerm;
    SOSLASSO.wholebrain.bootstrappedPerm};
% differences between methods
plot.meanDifference = [(mean(logLASSO.wholebrain.final) - mean(SOSLASSO.wholebrain.final));
    (mean(SOSLASSO.wholebrain.final) - mean(logLASSO.wholebrain.final))];
% bootstrapped permutation distribution for testing difference from the
% other method. N.B. These bootstrapped distributions are calculated as
% (neurally-inspired penalty - vanilla LASSO). Therefore, a positive
% difference score indicates better performance by the neurally-inspired
% penalty. To get the bootstrapped permutation distribution for testing
% whether vanilla LASSO is better, multiply this by -1. 
plot.bootstrappedPermDifference = {(classificationDifference.wholebrain.bootstrappedPerm)*-1;
    classificationDifference.wholebrain.bootstrappedPerm};
% do stats. For each row and column 
for r = 1:size(plot.mean,1)
    for c = 1:size(plot.mean,2)
        % difference from zero:
        % find the number of values in the bootstrapped permutation
        % distribution greater than the true mean
        b = sum(plot.bootstrappedPerm{r,c} > plot.mean(r,c));
        % calculate p-value using formula from Cox et al. (2024)
        plot.significance(r,c) = (b + 1)/10001;
        % difference between methods
        b = sum(plot.bootstrappedPermDifference{r,c} > plot.meanDifference(r,c));
        plot.significantDifference(r,c) = (b + 1)/10001;
    end
    % control the false discovery rate (Benjamini & Hochberg, 1995)
    plot.significance(r,:) = mafdr(plot.significance(r,:),'BHFDR',1);
    plot.significantDifference(r,:) = mafdr(plot.significantDifference(r,:),'BHFDR',1);
end

% do plots

% LASSO
analysisIndex = 1;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar(plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_wholebrain(1,:);
% set axes
xticks([])
yticks([])
% ylabel('Accuracy')
set(gca,'FontSize',20)
ylim([0.4,1])
% plot line at chance
yline(0.5,'--','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar(1:size(plot.mean,2),plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.025;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text((1:size(plot.mean,2))-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/logLASSO_wholebrain.svg','-dsvg')
close(fig)

% SOSLASSO 
analysisIndex = 2;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar(plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_wholebrain(1,:);
% set axes
xticks([])
yticks([])
% ylabel('Accuracy')
set(gca,'FontSize',20)
ylim([0.4,1])
% plot line at chance
yline(0.5,'--','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar(1:size(plot.mean,2),plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.025;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text((1:size(plot.mean,2))-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/SOSLASSO_wholebrain.svg','-dsvg')
close(fig)

% whole-brain RSL
clear plot
% prepare the variables - mean over participants
plot.mean = [mean(linearLASSO.wholebrain.D1.all.final), mean(linearLASSO.wholebrain.D2.all.final), mean(linearLASSO.wholebrain.D3.all.final); ...
    mean(grOWL.wholebrain.D1.all.final), mean(grOWL.wholebrain.D2.all.final), mean(grOWL.wholebrain.D3.all.final)];
% standard error
plot.standardError = [(std(linearLASSO.wholebrain.D1.all.final)/sqrt(size(linearLASSO.wholebrain.D1.all.final,2))), (std(linearLASSO.wholebrain.D2.all.final)/sqrt(size(linearLASSO.wholebrain.D2.all.final,2))), (std(linearLASSO.wholebrain.D3.all.final)/sqrt(size(linearLASSO.wholebrain.D3.all.final,2))); ...
    (std(grOWL.wholebrain.D1.all.final)/sqrt(size(grOWL.wholebrain.D1.all.final,2))), (std(grOWL.wholebrain.D2.all.final)/sqrt(size(grOWL.wholebrain.D2.all.final,2))), (std(grOWL.wholebrain.D3.all.final)/sqrt(size(grOWL.wholebrain.D3.all.final,2)))]; 
% bootstrapped permutation distribution for testing difference from zero
plot.bootstrappedPerm = {linearLASSO.wholebrain.D1.all.bootstrappedPerm, linearLASSO.wholebrain.D2.all.bootstrappedPerm,linearLASSO.wholebrain.D3.all.bootstrappedPerm; 
    grOWL.wholebrain.D1.all.bootstrappedPerm, grOWL.wholebrain.D2.all.bootstrappedPerm,grOWL.wholebrain.D3.all.bootstrappedPerm};
% differences between methods
plot.meanDifference = [(mean(linearLASSO.wholebrain.D1.all.final) - mean(grOWL.wholebrain.D1.all.final)), (mean(linearLASSO.wholebrain.D2.all.final) - mean(grOWL.wholebrain.D2.all.final)), (mean(linearLASSO.wholebrain.D3.all.final) - mean(grOWL.wholebrain.D3.all.final)); ...
    (mean(grOWL.wholebrain.D1.all.final) - mean(linearLASSO.wholebrain.D1.all.final)), (mean(grOWL.wholebrain.D2.all.final) - mean(linearLASSO.wholebrain.D2.all.final)), (mean(grOWL.wholebrain.D3.all.final) - mean(linearLASSO.wholebrain.D3.all.final))];
% bootstrapped permutation distribution for testing difference from the
% other method
plot.bootstrappedPermDifference = {(RSLDifference.wholebrain.D1.all.bootstrappedPerm)*-1, (RSLDifference.wholebrain.D2.all.bootstrappedPerm)*-1, (RSLDifference.wholebrain.D3.all.bootstrappedPerm)*-1; ...
    RSLDifference.wholebrain.D1.all.bootstrappedPerm, RSLDifference.wholebrain.D2.all.bootstrappedPerm, RSLDifference.wholebrain.D3.all.bootstrappedPerm};
% do stats. For each row and column 
for r = 1:size(plot.mean,1)
    for c = 1:size(plot.mean,2)
        % difference from zero:
        % find the number of values in the bootstrapped permutation
        % distribution greater than the true mean
        b = sum(plot.bootstrappedPerm{r,c} > plot.mean(r,c));
        % calculate p-value using formula from Cox et al. (2024)
        plot.significance(r,c) = (b + 1)/10001;
        % difference between methods
        b = sum(plot.bootstrappedPermDifference{r,c} > plot.meanDifference(r,c));
        plot.significantDifference(r,c) = (b + 1)/10001;
    end
    % control the false discovery rate (Benjamini & Hochberg, 1995)
    plot.significance(r,:) = mafdr(plot.significance(r,:),'BHFDR',1);
    plot.significantDifference(r,:) = mafdr(plot.significantDifference(r,:),'BHFDR',1);
end

% do plots

% LASSO
analysisIndex = 1;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar(plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_wholebrain(1,:);
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar(1:size(plot.mean,2),plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.05;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text((1:size(plot.mean,2))-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/linearLASSO_wholebrain_all.svg','-dsvg')
close(fig)

% grOWL
analysisIndex = 2;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar(plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_wholebrain(1,:);
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar(1:size(plot.mean,2),plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.05;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text((1:size(plot.mean,2))-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/grOWL_wholebrain_all.svg','-dsvg')
close(fig)

% animate only
clear plot
% prepare the variables - mean over participants
plot.mean = [mean(linearLASSO.wholebrain.D1.animate.final), mean(linearLASSO.wholebrain.D2.animate.final), mean(linearLASSO.wholebrain.D3.animate.final); ...
    mean(grOWL.wholebrain.D1.animate.final), mean(grOWL.wholebrain.D2.animate.final), mean(grOWL.wholebrain.D3.animate.final)];
% standard error
plot.standardError = [(std(linearLASSO.wholebrain.D1.animate.final)/sqrt(size(linearLASSO.wholebrain.D1.animate.final,2))), (std(linearLASSO.wholebrain.D2.animate.final)/sqrt(size(linearLASSO.wholebrain.D2.animate.final,2))), (std(linearLASSO.wholebrain.D3.animate.final)/sqrt(size(linearLASSO.wholebrain.D3.animate.final,2))); ...
    (std(grOWL.wholebrain.D1.animate.final)/sqrt(size(grOWL.wholebrain.D1.animate.final,2))), (std(grOWL.wholebrain.D2.animate.final)/sqrt(size(grOWL.wholebrain.D2.animate.final,2))), (std(grOWL.wholebrain.D3.animate.final)/sqrt(size(grOWL.wholebrain.D3.animate.final,2)))]; 
% bootstrapped permutation distribution for testing difference from zero
plot.bootstrappedPerm = {linearLASSO.wholebrain.D1.animate.bootstrappedPerm, linearLASSO.wholebrain.D2.animate.bootstrappedPerm,linearLASSO.wholebrain.D3.animate.bootstrappedPerm; 
    grOWL.wholebrain.D1.animate.bootstrappedPerm, grOWL.wholebrain.D2.animate.bootstrappedPerm,grOWL.wholebrain.D3.animate.bootstrappedPerm};
% differences between methods
plot.meanDifference = [(mean(linearLASSO.wholebrain.D1.animate.final) - mean(grOWL.wholebrain.D1.animate.final)), (mean(linearLASSO.wholebrain.D2.animate.final) - mean(grOWL.wholebrain.D2.animate.final)), (mean(linearLASSO.wholebrain.D3.animate.final) - mean(grOWL.wholebrain.D3.animate.final)); ...
    (mean(grOWL.wholebrain.D1.animate.final) - mean(linearLASSO.wholebrain.D1.animate.final)), (mean(grOWL.wholebrain.D2.animate.final) - mean(linearLASSO.wholebrain.D2.animate.final)), (mean(grOWL.wholebrain.D3.animate.final) - mean(linearLASSO.wholebrain.D3.animate.final))];
% bootstrapped permutation distribution for testing difference from the
% other method
plot.bootstrappedPermDifference = {(RSLDifference.wholebrain.D1.animate.bootstrappedPerm)*-1, (RSLDifference.wholebrain.D2.animate.bootstrappedPerm)*-1, (RSLDifference.wholebrain.D3.animate.bootstrappedPerm)*-1; ...
    RSLDifference.wholebrain.D1.animate.bootstrappedPerm, RSLDifference.wholebrain.D2.animate.bootstrappedPerm, RSLDifference.wholebrain.D3.animate.bootstrappedPerm};
% do stats. For each row and column 
for r = 1:size(plot.mean,1)
    for c = 1:size(plot.mean,2)
        % difference from zero:
        % find the number of values in the bootstrapped permutation
        % distribution greater than the true mean
        b = sum(plot.bootstrappedPerm{r,c} > plot.mean(r,c));
        % calculate p-value using formula from Cox et al. (2024)
        plot.significance(r,c) = (b + 1)/10001;
        % difference between methods
        b = sum(plot.bootstrappedPermDifference{r,c} > plot.meanDifference(r,c));
        plot.significantDifference(r,c) = (b + 1)/10001;
    end
    % control the false discovery rate (Benjamini & Hochberg, 1995)
    plot.significance(r,:) = mafdr(plot.significance(r,:),'BHFDR',1);
    plot.significantDifference(r,:) = mafdr(plot.significantDifference(r,:),'BHFDR',1);
end

% do plots

% LASSO
analysisIndex = 1;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar(plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_wholebrain(1,:);
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar(1:size(plot.mean,2),plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.05;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text((1:size(plot.mean,2))-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/linearLASSO_wholebrain_animate.svg','-dsvg')
close(fig)

% grOWL
analysisIndex = 2;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar(plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_wholebrain(1,:);
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar(1:size(plot.mean,2),plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.05;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text((1:size(plot.mean,2))-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/grOWL_wholebrain_animate.svg','-dsvg')
close(fig)

% inanimate only
clear plot
% prepare the variables - mean over participants
plot.mean = [mean(linearLASSO.wholebrain.D1.inanimate.final), mean(linearLASSO.wholebrain.D2.inanimate.final), mean(linearLASSO.wholebrain.D3.inanimate.final); ...
    mean(grOWL.wholebrain.D1.inanimate.final), mean(grOWL.wholebrain.D2.inanimate.final), mean(grOWL.wholebrain.D3.inanimate.final)];
% standard error
plot.standardError = [(std(linearLASSO.wholebrain.D1.inanimate.final)/sqrt(size(linearLASSO.wholebrain.D1.inanimate.final,2))), (std(linearLASSO.wholebrain.D2.inanimate.final)/sqrt(size(linearLASSO.wholebrain.D2.inanimate.final,2))), (std(linearLASSO.wholebrain.D3.inanimate.final)/sqrt(size(linearLASSO.wholebrain.D3.inanimate.final,2))); ...
    (std(grOWL.wholebrain.D1.inanimate.final)/sqrt(size(grOWL.wholebrain.D1.inanimate.final,2))), (std(grOWL.wholebrain.D2.inanimate.final)/sqrt(size(grOWL.wholebrain.D2.inanimate.final,2))), (std(grOWL.wholebrain.D3.inanimate.final)/sqrt(size(grOWL.wholebrain.D3.inanimate.final,2)))]; 
% bootstrapped permutation distribution for testing difference from zero
plot.bootstrappedPerm = {linearLASSO.wholebrain.D1.inanimate.bootstrappedPerm, linearLASSO.wholebrain.D2.inanimate.bootstrappedPerm,linearLASSO.wholebrain.D3.inanimate.bootstrappedPerm; 
    grOWL.wholebrain.D1.inanimate.bootstrappedPerm, grOWL.wholebrain.D2.inanimate.bootstrappedPerm,grOWL.wholebrain.D3.inanimate.bootstrappedPerm};
% differences between methods
plot.meanDifference = [(mean(linearLASSO.wholebrain.D1.inanimate.final) - mean(grOWL.wholebrain.D1.inanimate.final)), (mean(linearLASSO.wholebrain.D2.inanimate.final) - mean(grOWL.wholebrain.D2.inanimate.final)), (mean(linearLASSO.wholebrain.D3.inanimate.final) - mean(grOWL.wholebrain.D3.inanimate.final)); ...
    (mean(grOWL.wholebrain.D1.inanimate.final) - mean(linearLASSO.wholebrain.D1.inanimate.final)), (mean(grOWL.wholebrain.D2.inanimate.final) - mean(linearLASSO.wholebrain.D2.inanimate.final)), (mean(grOWL.wholebrain.D3.inanimate.final) - mean(linearLASSO.wholebrain.D3.inanimate.final))];
% bootstrapped permutation distribution for testing difference from the
% other method
plot.bootstrappedPermDifference = {(RSLDifference.wholebrain.D1.inanimate.bootstrappedPerm)*-1, (RSLDifference.wholebrain.D2.inanimate.bootstrappedPerm)*-1, (RSLDifference.wholebrain.D3.inanimate.bootstrappedPerm)*-1; ...
    RSLDifference.wholebrain.D1.inanimate.bootstrappedPerm, RSLDifference.wholebrain.D2.inanimate.bootstrappedPerm, RSLDifference.wholebrain.D3.inanimate.bootstrappedPerm};
% do stats. For each row and column 
for r = 1:size(plot.mean,1)
    for c = 1:size(plot.mean,2)
        % difference from zero:
        % find the number of values in the bootstrapped permutation
        % distribution greater than the true mean
        b = sum(plot.bootstrappedPerm{r,c} > plot.mean(r,c));
        % calculate p-value using formula from Cox et al. (2024)
        plot.significance(r,c) = (b + 1)/10001;
        % difference between methods
        b = sum(plot.bootstrappedPermDifference{r,c} > plot.meanDifference(r,c));
        plot.significantDifference(r,c) = (b + 1)/10001;
    end
    % control the false discovery rate (Benjamini & Hochberg, 1995)
    plot.significance(r,:) = mafdr(plot.significance(r,:),'BHFDR',1);
    plot.significantDifference(r,:) = mafdr(plot.significantDifference(r,:),'BHFDR',1);
end

% do plots

% LASSO
analysisIndex = 1;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar(plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_wholebrain(1,:);
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar(1:size(plot.mean,2),plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.05;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text((1:size(plot.mean,2))-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/linearLASSO_wholebrain_inanimate.svg','-dsvg')
close(fig)

% grOWL
analysisIndex = 2;
fig = figure;
% make it large
fig.Position = [610 240 900 725];
% plot bars
b = bar(plot.mean(analysisIndex,:));
% set colours
b.FaceColor = 'flat';
b.EdgeColor = 'none';
b.CData = colours.colours_LASSO_wholebrain(1,:);
% set axes
xticks([])
yticks([])
% set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
% ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
% plot thicker x-axis
yline(0,'-','LineWidth',3,'Alpha',1)
% plot standard error
hold on
errorbar(1:size(plot.mean,2),plot.mean(analysisIndex,:),plot.standardError(analysisIndex,:),'LineStyle','none','LineWidth',3,'CapSize',30,'Color','black');
% make significance stars
stars = {};
for c = 1:size(plot.significance,2)
    if plot.significance(analysisIndex,c) < 0.05 && plot.significantDifference(analysisIndex,c) < 0.05
        stars{c} = '*+';
    elseif plot.significance(analysisIndex,c) < 0.05
        stars{c} = '*';
    elseif plot.significance(analysisIndex,c) >= 0.05
       stars{c} = '';
    end
end
% calculate where the stars should go on the page
positionOnPage = plot.mean(analysisIndex,:) + plot.standardError(analysisIndex,:) + 0.05;
positionOnPage(positionOnPage > 0.975) = 0.975;
% plot the stars
text((1:size(plot.mean,2))-0.1,positionOnPage,stars,'FontSize',40,'FontWeight','Bold');
% save figure in a way that maintains its high resolution
set(fig,'Renderer','painters')
print(fig,'/group/mlr-lab/Saskia/7T_decoding/20260122_bar_charts/grOWL_wholebrain_inanimate.svg','-dsvg')
close(fig)

% finally, conduct a statistical test for dynamism. 
clear plot
% prepare the variables - differences between regions for each method
plot.meanDifference = [mean(dynamismDifference.LASSO.left.final); ...
    mean(dynamismDifference.LASSO.right.final); ...
    mean(dynamismDifference.SOSLASSO.left.final); ...
    mean(dynamismDifference.SOSLASSO.right.final)];
% bootstrapped permutation distribution for testing differences between
% regions
plot.bootstrappedPermDifference = {dynamismDifference.LASSO.left.bootstrappedPerm; ...
    dynamismDifference.LASSO.right.bootstrappedPerm; ...
    dynamismDifference.SOSLASSO.left.bootstrappedPerm; ...
    dynamismDifference.SOSLASSO.right.bootstrappedPerm};
% do stats. for each row
for r = 1:size(plot.meanDifference,1)
    % differences between methods
    b = sum(plot.bootstrappedPermDifference{r,1} > plot.meanDifference(r,1));
    plot.significantDifference(r,1) = (b + 1)/10001;
end

