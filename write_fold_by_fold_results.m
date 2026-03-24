% collates results and writes correlations for each penalty, participant,
% and fold. For troubleshooting statistics!

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

    % organise permutation results: seed x participant x fold
    logLASSO.wholebrain.perm(:,s,:) = tmp.averaged_whole_brain.permoutput';
    logLASSO.LATL.perm(:,s,:) = tmp.averaged_left_ATL.permoutput';
    logLASSO.LATLant.perm(:,s,:) = tmp.averaged_left_ATL_ant.permoutput';
    logLASSO.LATLpos.perm(:,s,:) = tmp.averaged_left_ATL_pos.permoutput';
    logLASSO.RATL.perm(:,s,:) = tmp.averaged_right_ATL.permoutput';
    logLASSO.RATLant.perm(:,s,:) = tmp.averaged_right_ATL_ant.permoutput';
    logLASSO.RATLpos.perm(:,s,:) = tmp.averaged_right_ATL_pos.permoutput';
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

% write .csvs. Final results first

% set up table in which to write classification results
ROINames = {'wholebrain', 'LATL', 'LATLant', 'LATLpos', 'RATL','RATLant', 'RATLpos'};

penalty = repmat('LASSO',7*size(subcode,2)*10,1); % 7 ROIs x 10 folds
ROI = repelem(ROINames',size(subcode,2)*10,1);
subject = repmat(repelem(subcode',10,1),7,1);
holdout = repmat((1:10)',size(subcode,2)*7,1);
resultsTable = table(penalty,ROI,subject,holdout);

% fill in LASSO results
resultsTable.AUC = [reshape(logLASSO.wholebrain.final,[],1); reshape(logLASSO.LATL.final,[],1); reshape(logLASSO.LATLant.final,[],1); reshape(logLASSO.LATLpos.final,[],1); reshape(logLASSO.RATL.final,[],1); reshape(logLASSO.RATLant.final,[],1); reshape(logLASSO.RATLpos.final,[],1)];
% save
if ~exist([fileparts(root),'/work/fold_by_fold_results/'])
    mkdir([fileparts(root),'/work/fold_by_fold_results/']);
end
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/classification_LASSO_final.csv'],'Delimiter',',');

% overwrite with SOSLASSO results
resultsTable.penalty = repmat('SOSLASSO',7*size(subcode,2)*10,1); % 7 ROIs x 10 folds
resultsTable.AUC = [reshape(SOSLASSO.wholebrain.final,[],1); reshape(SOSLASSO.LATL.final,[],1); reshape(SOSLASSO.LATLant.final,[],1); reshape(SOSLASSO.LATLpos.final,[],1); reshape(SOSLASSO.RATL.final,[],1); reshape(SOSLASSO.RATLant.final,[],1); reshape(SOSLASSO.RATLpos.final,[],1)];
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/classification_SOSLASSO_final.csv'],'Delimiter',',');

% set up table in which to write RSL results
clear resultsTable
penalty = repmat('LASSO',7*size(subcode,2)*3*10,1); % 7 ROIs x 3 dimensions x 10 folds
items = repmat('all',7*size(subcode,2)*3*10,1);
ROI = repelem(ROINames',size(subcode,2)*3*10,1);
dimension = repmat(repelem((1:3)',size(subcode,2)*10,1),7,1);
subject = repmat(repelem(subcode',10,1),3*7,1);
holdout = repmat((1:10)',size(subcode,2)*3*7,1);
resultsTable = table(penalty,items,ROI,dimension,subject,holdout);

% fill in LASSO results for all stimuli
resultsTable.correlation = [reshape(linearLASSO.wholebrain.D1.all.final,[],1); reshape(linearLASSO.wholebrain.D2.all.final,[],1); reshape(linearLASSO.wholebrain.D3.all.final,[],1); ...
    reshape(linearLASSO.LATL.D1.all.final,[],1); reshape(linearLASSO.LATL.D2.all.final,[],1); reshape(linearLASSO.LATL.D3.all.final,[],1); ...
    reshape(linearLASSO.LATLant.D1.all.final,[],1); reshape(linearLASSO.LATLant.D2.all.final,[],1); reshape(linearLASSO.LATLant.D3.all.final,[],1); ...
    reshape(linearLASSO.LATLpos.D1.all.final,[],1); reshape(linearLASSO.LATLpos.D2.all.final,[],1); reshape(linearLASSO.LATLpos.D3.all.final,[],1); ...
    reshape(linearLASSO.RATL.D1.all.final,[],1); reshape(linearLASSO.RATL.D2.all.final,[],1); reshape(linearLASSO.RATL.D3.all.final,[],1); ...
    reshape(linearLASSO.RATLant.D1.all.final,[],1); reshape(linearLASSO.RATLant.D2.all.final,[],1); reshape(linearLASSO.RATLant.D3.all.final,[],1); ...
    reshape(linearLASSO.RATLpos.D1.all.final,[],1); reshape(linearLASSO.RATLpos.D2.all.final,[],1); reshape(linearLASSO.RATLpos.D3.all.final,[],1)]; 
% save
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/RSL_LASSO_all_final.csv'],'Delimiter',',');

% overwrite with results for animate stimuli
resultsTable.items = repmat('animate',7*size(subcode,2)*3*10,1);
resultsTable.correlation = [reshape(linearLASSO.wholebrain.D1.animate.final,[],1); reshape(linearLASSO.wholebrain.D2.animate.final,[],1); reshape(linearLASSO.wholebrain.D3.animate.final,[],1); ...
    reshape(linearLASSO.LATL.D1.animate.final,[],1); reshape(linearLASSO.LATL.D2.animate.final,[],1); reshape(linearLASSO.LATL.D3.animate.final,[],1); ...
    reshape(linearLASSO.LATLant.D1.animate.final,[],1); reshape(linearLASSO.LATLant.D2.animate.final,[],1); reshape(linearLASSO.LATLant.D3.animate.final,[],1); ...
    reshape(linearLASSO.LATLpos.D1.animate.final,[],1); reshape(linearLASSO.LATLpos.D2.animate.final,[],1); reshape(linearLASSO.LATLpos.D3.animate.final,[],1); ...
    reshape(linearLASSO.RATL.D1.animate.final,[],1); reshape(linearLASSO.RATL.D2.animate.final,[],1); reshape(linearLASSO.RATL.D3.animate.final,[],1); ...
    reshape(linearLASSO.RATLant.D1.animate.final,[],1); reshape(linearLASSO.RATLant.D2.animate.final,[],1); reshape(linearLASSO.RATLant.D3.animate.final,[],1); ...
    reshape(linearLASSO.RATLpos.D1.animate.final,[],1); reshape(linearLASSO.RATLpos.D2.animate.final,[],1); reshape(linearLASSO.RATLpos.D3.animate.final,[],1)]; 
% save
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/RSL_LASSO_animate_final.csv'],'Delimiter',',');

% overwrite with results for inanimate stimuli
resultsTable.items = repmat('inanimate',7*size(subcode,2)*3*10,1);
resultsTable.correlation = [reshape(linearLASSO.wholebrain.D1.inanimate.final,[],1); reshape(linearLASSO.wholebrain.D2.inanimate.final,[],1); reshape(linearLASSO.wholebrain.D3.inanimate.final,[],1); ...
    reshape(linearLASSO.LATL.D1.inanimate.final,[],1); reshape(linearLASSO.LATL.D2.inanimate.final,[],1); reshape(linearLASSO.LATL.D3.inanimate.final,[],1); ...
    reshape(linearLASSO.LATLant.D1.inanimate.final,[],1); reshape(linearLASSO.LATLant.D2.inanimate.final,[],1); reshape(linearLASSO.LATLant.D3.inanimate.final,[],1); ...
    reshape(linearLASSO.LATLpos.D1.inanimate.final,[],1); reshape(linearLASSO.LATLpos.D2.inanimate.final,[],1); reshape(linearLASSO.LATLpos.D3.inanimate.final,[],1); ...
    reshape(linearLASSO.RATL.D1.inanimate.final,[],1); reshape(linearLASSO.RATL.D2.inanimate.final,[],1); reshape(linearLASSO.RATL.D3.inanimate.final,[],1); ...
    reshape(linearLASSO.RATLant.D1.inanimate.final,[],1); reshape(linearLASSO.RATLant.D2.inanimate.final,[],1); reshape(linearLASSO.RATLant.D3.inanimate.final,[],1); ...
    reshape(linearLASSO.RATLpos.D1.inanimate.final,[],1); reshape(linearLASSO.RATLpos.D2.inanimate.final,[],1); reshape(linearLASSO.RATLpos.D3.inanimate.final,[],1)]; 
% save
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/RSL_LASSO_inanimate_final.csv'],'Delimiter',',');

% overwrite with grOWL results for all stimuli
resultsTable.penalty = repmat('grOWL',7*size(subcode,2)*3*10,1);
resultsTable.items = repmat('all',7*size(subcode,2)*3*10,1);
resultsTable.correlation = [reshape(grOWL.wholebrain.D1.all.final,[],1); reshape(grOWL.wholebrain.D2.all.final,[],1); reshape(grOWL.wholebrain.D3.all.final,[],1); ...
    reshape(grOWL.LATL.D1.all.final,[],1); reshape(grOWL.LATL.D2.all.final,[],1); reshape(grOWL.LATL.D3.all.final,[],1); ...
    reshape(grOWL.LATLant.D1.all.final,[],1); reshape(grOWL.LATLant.D2.all.final,[],1); reshape(grOWL.LATLant.D3.all.final,[],1); ...
    reshape(grOWL.LATLpos.D1.all.final,[],1); reshape(grOWL.LATLpos.D2.all.final,[],1); reshape(grOWL.LATLpos.D3.all.final,[],1); ...
    reshape(grOWL.RATL.D1.all.final,[],1); reshape(grOWL.RATL.D2.all.final,[],1); reshape(grOWL.RATL.D3.all.final,[],1); ...
    reshape(grOWL.RATLant.D1.all.final,[],1); reshape(grOWL.RATLant.D2.all.final,[],1); reshape(grOWL.RATLant.D3.all.final,[],1); ...
    reshape(grOWL.RATLpos.D1.all.final,[],1); reshape(grOWL.RATLpos.D2.all.final,[],1); reshape(grOWL.RATLpos.D3.all.final,[],1)]; 
% save
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/RSL_grOWL_all_final.csv'],'Delimiter',',');

% overwrite with results for animate stimuli
resultsTable.items = repmat('animate',7*size(subcode,2)*3*10,1);
resultsTable.correlation = [reshape(grOWL.wholebrain.D1.animate.final,[],1); reshape(grOWL.wholebrain.D2.animate.final,[],1); reshape(grOWL.wholebrain.D3.animate.final,[],1); ...
    reshape(grOWL.LATL.D1.animate.final,[],1); reshape(grOWL.LATL.D2.animate.final,[],1); reshape(grOWL.LATL.D3.animate.final,[],1); ...
    reshape(grOWL.LATLant.D1.animate.final,[],1); reshape(grOWL.LATLant.D2.animate.final,[],1); reshape(grOWL.LATLant.D3.animate.final,[],1); ...
    reshape(grOWL.LATLpos.D1.animate.final,[],1); reshape(grOWL.LATLpos.D2.animate.final,[],1); reshape(grOWL.LATLpos.D3.animate.final,[],1); ...
    reshape(grOWL.RATL.D1.animate.final,[],1); reshape(grOWL.RATL.D2.animate.final,[],1); reshape(grOWL.RATL.D3.animate.final,[],1); ...
    reshape(grOWL.RATLant.D1.animate.final,[],1); reshape(grOWL.RATLant.D2.animate.final,[],1); reshape(grOWL.RATLant.D3.animate.final,[],1); ...
    reshape(grOWL.RATLpos.D1.animate.final,[],1); reshape(grOWL.RATLpos.D2.animate.final,[],1); reshape(grOWL.RATLpos.D3.animate.final,[],1)]; 
% save
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/RSL_grOWL_animate_final.csv'],'Delimiter',',');

% overwrite with results for inanimate stimuli
resultsTable.items = repmat('inanimate',7*size(subcode,2)*3*10,1);
resultsTable.correlation = [reshape(grOWL.wholebrain.D1.inanimate.final,[],1); reshape(grOWL.wholebrain.D2.inanimate.final,[],1); reshape(grOWL.wholebrain.D3.inanimate.final,[],1); ...
    reshape(grOWL.LATL.D1.inanimate.final,[],1); reshape(grOWL.LATL.D2.inanimate.final,[],1); reshape(grOWL.LATL.D3.inanimate.final,[],1); ...
    reshape(grOWL.LATLant.D1.inanimate.final,[],1); reshape(grOWL.LATLant.D2.inanimate.final,[],1); reshape(grOWL.LATLant.D3.inanimate.final,[],1); ...
    reshape(grOWL.LATLpos.D1.inanimate.final,[],1); reshape(grOWL.LATLpos.D2.inanimate.final,[],1); reshape(grOWL.LATLpos.D3.inanimate.final,[],1); ...
    reshape(grOWL.RATL.D1.inanimate.final,[],1); reshape(grOWL.RATL.D2.inanimate.final,[],1); reshape(grOWL.RATL.D3.inanimate.final,[],1); ...
    reshape(grOWL.RATLant.D1.inanimate.final,[],1); reshape(grOWL.RATLant.D2.inanimate.final,[],1); reshape(grOWL.RATLant.D3.inanimate.final,[],1); ...
    reshape(grOWL.RATLpos.D1.inanimate.final,[],1); reshape(grOWL.RATLpos.D2.inanimate.final,[],1); reshape(grOWL.RATLpos.D3.inanimate.final,[],1)]; 
% save
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/RSL_grOWL_inanimate_final.csv'],'Delimiter',',');

% now permutation results

% set up table in which to write classification results
clear resultsTable
penalty = repmat('LASSO',7*size(subcode,2)*10*100,1); % 7 ROIs x 10 folds x 100 permutations
ROI = repelem(ROINames',size(subcode,2)*10*100,1); 
subject = repmat(repelem(subcode',10*100,1),7,1);
holdout = repmat(repelem((1:10)',100),size(subcode,2)*7,1);
randomseed = repmat((1:100)',7*10*size(subcode,2),1);
resultsTable = table(penalty,ROI,subject,holdout,randomseed);

% fill in LASSO results
resultsTable.AUC = [reshape(permute(logLASSO.wholebrain.perm,[1,3,2]),[],1); reshape(permute(logLASSO.LATL.perm,[1,3,2]),[],1); reshape(permute(logLASSO.LATLant.perm,[1,3,2]),[],1); reshape(permute(logLASSO.LATLpos.perm,[1,3,2]),[],1); reshape(permute(logLASSO.RATL.perm,[1,3,2]),[],1); reshape(permute(logLASSO.RATLant.perm,[1,3,2]),[],1); reshape(permute(logLASSO.RATLpos.perm,[1,3,2]),[],1)];
% save
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/classification_LASSO_perm.csv'],'Delimiter',',');

% overwrite with SOSLASSO results
resultsTable.penalty = repmat('SOSLASSO',7*size(subcode,2)*10*100,1);
resultsTable.AUC = [reshape(permute(SOSLASSO.wholebrain.perm,[1,3,2]),[],1); reshape(permute(SOSLASSO.LATL.perm,[1,3,2]),[],1); reshape(permute(SOSLASSO.LATLant.perm,[1,3,2]),[],1); reshape(permute(SOSLASSO.LATLpos.perm,[1,3,2]),[],1); reshape(permute(SOSLASSO.RATL.perm,[1,3,2]),[],1); reshape(permute(SOSLASSO.RATLant.perm,[1,3,2]),[],1); reshape(permute(SOSLASSO.RATLpos.perm,[1,3,2]),[],1)];
% save
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/classification_SOSLASSO_perm.csv'],'Delimiter',',');

% set up table in which to write RSL results
clear resultsTable
penalty = repmat('LASSO',7*size(subcode,2)*3*10*100,1); % 7 ROIs x 3 dimensions x 10 folds x 100 permutations
items = repmat('all',7*size(subcode,2)*3*10*100,1);
ROI = repelem(ROINames',size(subcode,2)*3*10*100,1);
dimension = repmat(repelem((1:3)',size(subcode,2)*10*100,1),7,1);
subject = repmat(repelem(subcode',10*100,1),3*7,1);
holdout = repmat(repelem((1:10)',100),size(subcode,2)*3*7,1);
randomseed = repmat((1:100)',7*3*10*size(subcode,2),1);
resultsTable = table(penalty,items,ROI,dimension,subject,holdout,randomseed);

% fill in LASSO results for all stimuli
resultsTable.correlation = [reshape(permute(linearLASSO.wholebrain.D1.all.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.wholebrain.D2.all.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.wholebrain.D3.all.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.LATL.D1.all.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATL.D2.all.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATL.D3.all.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.LATLant.D1.all.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATLant.D2.all.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATLant.D3.all.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.LATLpos.D1.all.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATLpos.D2.all.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATLpos.D3.all.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.RATL.D1.all.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATL.D2.all.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATL.D3.all.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.RATLant.D1.all.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATLant.D2.all.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATLant.D3.all.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.RATLpos.D1.all.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATLpos.D2.all.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATLpos.D3.all.perm,[1,3,2]),[],1)]; 
% save
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/RSL_LASSO_all_perm.csv'],'Delimiter',',');

% overwrite with results for animate stimuli
resultsTable.items = repmat('animate',7*size(subcode,2)*3*10*100,1);
resultsTable.correlation = [reshape(permute(linearLASSO.wholebrain.D1.animate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.wholebrain.D2.animate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.wholebrain.D3.animate.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.LATL.D1.animate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATL.D2.animate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATL.D3.animate.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.LATLant.D1.animate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATLant.D2.animate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATLant.D3.animate.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.LATLpos.D1.animate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATLpos.D2.animate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATLpos.D3.animate.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.RATL.D1.animate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATL.D2.animate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATL.D3.animate.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.RATLant.D1.animate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATLant.D2.animate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATLant.D3.animate.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.RATLpos.D1.animate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATLpos.D2.animate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATLpos.D3.animate.perm,[1,3,2]),[],1)]; 
% save
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/RSL_LASSO_animate_perm.csv'],'Delimiter',',');

% overwrite with results for inanimate stimuli
resultsTable.items = repmat('inanimate',7*size(subcode,2)*3*10*100,1);
resultsTable.correlation = [reshape(permute(linearLASSO.wholebrain.D1.inanimate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.wholebrain.D2.inanimate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.wholebrain.D3.inanimate.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.LATL.D1.inanimate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATL.D2.inanimate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATL.D3.inanimate.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.LATLant.D1.inanimate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATLant.D2.inanimate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATLant.D3.inanimate.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.LATLpos.D1.inanimate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATLpos.D2.inanimate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.LATLpos.D3.inanimate.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.RATL.D1.inanimate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATL.D2.inanimate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATL.D3.inanimate.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.RATLant.D1.inanimate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATLant.D2.inanimate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATLant.D3.inanimate.perm,[1,3,2]),[],1); ...
    reshape(permute(linearLASSO.RATLpos.D1.inanimate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATLpos.D2.inanimate.perm,[1,3,2]),[],1); reshape(permute(linearLASSO.RATLpos.D3.inanimate.perm,[1,3,2]),[],1)]; 
% save
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/RSL_LASSO_inanimate_perm.csv'],'Delimiter',',');

% overwrite with grOWL results for all stimuli
resultsTable.penalty = repmat('grOWL',7*size(subcode,2)*3*10*100,1);
resultsTable.items = repmat('all',7*size(subcode,2)*3*10*100,1);
resultsTable.correlation = [reshape(permute(grOWL.wholebrain.D1.all.perm,[1,3,2]),[],1); reshape(permute(grOWL.wholebrain.D2.all.perm,[1,3,2]),[],1); reshape(permute(grOWL.wholebrain.D3.all.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.LATL.D1.all.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATL.D2.all.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATL.D3.all.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.LATLant.D1.all.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATLant.D2.all.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATLant.D3.all.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.LATLpos.D1.all.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATLpos.D2.all.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATLpos.D3.all.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.RATL.D1.all.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATL.D2.all.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATL.D3.all.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.RATLant.D1.all.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATLant.D2.all.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATLant.D3.all.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.RATLpos.D1.all.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATLpos.D2.all.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATLpos.D3.all.perm,[1,3,2]),[],1)]; 
% save
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/RSL_grOWL_all_perm.csv'],'Delimiter',',');

% overwrite with results for animate stimuli
resultsTable.items = repmat('animate',7*size(subcode,2)*3*10*100,1);
resultsTable.correlation = [reshape(permute(grOWL.wholebrain.D1.animate.perm,[1,3,2]),[],1); reshape(permute(grOWL.wholebrain.D2.animate.perm,[1,3,2]),[],1); reshape(permute(grOWL.wholebrain.D3.animate.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.LATL.D1.animate.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATL.D2.animate.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATL.D3.animate.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.LATLant.D1.animate.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATLant.D2.animate.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATLant.D3.animate.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.LATLpos.D1.animate.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATLpos.D2.animate.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATLpos.D3.animate.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.RATL.D1.animate.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATL.D2.animate.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATL.D3.animate.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.RATLant.D1.animate.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATLant.D2.animate.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATLant.D3.animate.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.RATLpos.D1.animate.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATLpos.D2.animate.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATLpos.D3.animate.perm,[1,3,2]),[],1)]; 
% save
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/RSL_grOWL_animate_perm.csv'],'Delimiter',',');

% overwrite with results for inanimate stimuli
resultsTable.items = repmat('inanimate',7*size(subcode,2)*3*10*100,1);
resultsTable.correlation = [reshape(permute(grOWL.wholebrain.D1.inanimate.perm,[1,3,2]),[],1); reshape(permute(grOWL.wholebrain.D2.inanimate.perm,[1,3,2]),[],1); reshape(permute(grOWL.wholebrain.D3.inanimate.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.LATL.D1.inanimate.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATL.D2.inanimate.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATL.D3.inanimate.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.LATLant.D1.inanimate.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATLant.D2.inanimate.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATLant.D3.inanimate.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.LATLpos.D1.inanimate.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATLpos.D2.inanimate.perm,[1,3,2]),[],1); reshape(permute(grOWL.LATLpos.D3.inanimate.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.RATL.D1.inanimate.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATL.D2.inanimate.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATL.D3.inanimate.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.RATLant.D1.inanimate.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATLant.D2.inanimate.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATLant.D3.inanimate.perm,[1,3,2]),[],1); ...
    reshape(permute(grOWL.RATLpos.D1.inanimate.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATLpos.D2.inanimate.perm,[1,3,2]),[],1); reshape(permute(grOWL.RATLpos.D3.inanimate.perm,[1,3,2]),[],1)]; 
% save
writetable(resultsTable,[fileparts(root),'/work/fold_by_fold_results/RSL_grOWL_inanimate_perm.csv'],'Delimiter',',');
