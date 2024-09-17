% by Saskia. Collates results and makes bar charts for thesis and
% manuscript writing. 

addpath('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/');
addpath('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/WISC_MVPA/src/');
root = ['/imaging/projects/cbu/wbic-p00591-DAISY/main/derivatives'];
cd([root]);

% get metadata
load([root,'/cox/tedana_averaged_metadata.mat'])

subcode = {'sub-001','sub-002','sub-003','sub-004','sub-007','sub-009','sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019','sub-020','sub-021','sub-022','sub-023','sub-024','sub-026','sub-028','sub-029','sub-030','sub-031','sub-032'};
% labelling of folders on condor
ROI = {'wholebrain','LATL','LATLant','LATLpos','RATL','RATLant','RATLpos'};

% define colour palettes
load(['/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/colour_palettes.mat'])

% 1. CLASSIFICATION (Logistic LASSO and SOSLASSO)

% 1A. Organise logistic LASSO results
whole_brain_log_LASSO = zeros(10,size(subcode,2));
left_ATL_log_LASSO = zeros(10,size(subcode,2));
right_ATL_log_LASSO = zeros(10,size(subcode,2));
left_ATL_ant_log_LASSO = zeros(10,size(subcode,2));
left_ATL_pos_log_LASSO = zeros(10,size(subcode,2));
right_ATL_ant_log_LASSO = zeros(10,size(subcode,2));
right_ATL_pos_log_LASSO = zeros(10,size(subcode,2));
% perm results are arranged !! participant x random seed
perm_whole_brain_log_LASSO = zeros(size(subcode,2),100);
perm_left_ATL_log_LASSO = zeros(size(subcode,2),100);
perm_right_ATL_log_LASSO = zeros(size(subcode,2),100);
perm_left_ATL_ant_log_LASSO = zeros(size(subcode,2),100);
perm_left_ATL_pos_log_LASSO = zeros(size(subcode,2),100);
perm_right_ATL_ant_log_LASSO = zeros(size(subcode,2),100);
perm_right_ATL_pos_log_LASSO = zeros(size(subcode,2),100);

% for every participant
for s = 1:size(subcode,2)

    % load results
    result = load([root,'/LASSO/log/tedana/',subcode{s},'/results.mat']);

    % collate 
    whole_brain_log_LASSO(:,s) = result.averaged_whole_brain.output;
    left_ATL_log_LASSO(:,s) = result.averaged_left_ATL.output;
    right_ATL_log_LASSO(:,s) = result.averaged_right_ATL.output;
    left_ATL_ant_log_LASSO(:,s) = result.averaged_left_ATL_ant.output;
    left_ATL_pos_log_LASSO(:,s) = result.averaged_left_ATL_pos.output;
    right_ATL_ant_log_LASSO(:,s) = result.averaged_right_ATL_ant.output;
    right_ATL_pos_log_LASSO(:,s) = result.averaged_right_ATL_pos.output;
 
    perm_whole_brain_log_LASSO(s,:) = mean(result.averaged_whole_brain.permoutput)';
    perm_left_ATL_log_LASSO(s,:) = mean(result.averaged_left_ATL.permoutput)';
    perm_right_ATL_log_LASSO(s,:) = mean(result.averaged_right_ATL.permoutput)';
    perm_left_ATL_ant_log_LASSO(s,:) = mean(result.averaged_left_ATL_ant.permoutput)';
    perm_left_ATL_pos_log_LASSO(s,:) = mean(result.averaged_left_ATL_pos.permoutput)';
    perm_right_ATL_ant_log_LASSO(s,:) = mean(result.averaged_right_ATL_ant.permoutput)';
    perm_right_ATL_pos_log_LASSO(s,:) = mean(result.averaged_right_ATL_pos.permoutput)';

end

% 1B. Organise SOSLASSO results

% initialise temporary results
result = {};
permresult = {};

% load results from condor

% for each ROI
for r = 1:length(ROI)
    
    % load results
    load([root,'/condor/derivatives/classification/SOSLASSO/performance/',ROI{r},'/final_performance.mat']);
    
    % initialise output
    output = zeros(10,size(subcode,2));

    % for each participant
    for s = 1:size(subcode,2)
        % get just the rows applying to that participant
        tmp = str2num(erase(subcode{s},'sub-'));
        subtable = Tallcv(Tallcv.subject==tmp,:);
        
        % calculate AUC for each fold.
        for ho = 1:10 
            % get target for that holdout set
            ytst = [zeros(5,1);ones(5,1)]; 
            % get predicted coordinates
            allpred = subtable.Yz{ho,1};
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate AUC
            [~,~,~,AUC] = perfcurve(ytst,pred,1);
            output(ho,s) = AUC;
        end
    end
    result(r).output = output;

    % load perm results
    load([root,'/condor/derivatives/classification/SOSLASSO/performance/',ROI{r},'/perm_performance.mat']);
    
    % initialise output - !! participant x random seed (x fold, but we will
    % average over folds)
    output = zeros(size(subcode,2), 100, 10);

    % for each participant
    for s = 1:size(subcode,2)
        % get just the rows applying to that participant
        tmp = str2num(erase(subcode{s},'sub-'));
        subtable = Tallcv(Tallcv.subject==tmp,:);
        
        % for every random seed
        for seed = 1:100
            subsubtable = subtable(subtable.RandomSeed == seed,:);
            % calculate AUC for each fold.
            for ho = 1:10 
                % get target for that holdout set
                ytst = [zeros(5,1);ones(5,1)]; 
                % get predicted coordinates
                allpred = subsubtable.Yz{ho,1};
                % index just coordinates in that holdout set
                pred = allpred(metadata(s).cvind(:,1)==ho);
                % calculate AUC
                [~,~,~,AUC] = perfcurve(ytst,pred,1);
                output(s,seed,ho) = AUC;
            end
        end
    end
    permresult(r).output = output;
end

% bring output into workspace
whole_brain_SOSLASSO = result(1).output;
left_ATL_SOSLASSO = result(2).output;
left_ATL_ant_SOSLASSO = result(3).output;
left_ATL_pos_SOSLASSO = result(4).output;
right_ATL_SOSLASSO = result(5).output;
right_ATL_ant_SOSLASSO = result(6).output;
right_ATL_pos_SOSLASSO = result(7).output;
perm_whole_brain_SOSLASSO = mean(permresult(1).output,3);
perm_left_ATL_SOSLASSO = mean(permresult(2).output,3);
perm_left_ATL_ant_SOSLASSO = mean(permresult(3).output,3);
perm_left_ATL_pos_SOSLASSO = mean(permresult(4).output,3);
perm_right_ATL_SOSLASSO = mean(permresult(5).output,3);
perm_right_ATL_ant_SOSLASSO = mean(permresult(6).output,3);
perm_right_ATL_pos_SOSLASSO = mean(permresult(7).output,3);

% 1C. Statistics with Stelzer bootstrapping for the ATL

% initialise collated results
group_log_LASSO_left = zeros(size(subcode,2),3);
group_log_LASSO_right = zeros(size(subcode,2),3);
group_SOSLASSO_left = zeros(size(subcode,2),3);
group_SOSLASSO_right = zeros(size(subcode,2),3);

% collate
group_log_LASSO_left(:,1) = mean(left_ATL_log_LASSO)';
group_log_LASSO_left(:,2) = mean(left_ATL_ant_log_LASSO)';
group_log_LASSO_left(:,3) = mean(left_ATL_pos_log_LASSO)';
group_log_LASSO_right(:,1) = mean(right_ATL_log_LASSO)';
group_log_LASSO_right(:,2) = mean(right_ATL_ant_log_LASSO)';
group_log_LASSO_right(:,3) = mean(right_ATL_pos_log_LASSO)';
group_SOSLASSO_left(:,1) = mean(left_ATL_SOSLASSO)';
group_SOSLASSO_left(:,2) = mean(left_ATL_ant_SOSLASSO)';
group_SOSLASSO_left(:,3) = mean(left_ATL_pos_SOSLASSO)';
group_SOSLASSO_right(:,1) = mean(right_ATL_SOSLASSO)';
group_SOSLASSO_right(:,2) = mean(right_ATL_ant_SOSLASSO)';
group_SOSLASSO_right(:,3) = mean(right_ATL_pos_SOSLASSO)';

% construct permutation distributions by randomly sampling group averages
% (see Cox et al., 2024, Imaging Neuroscience for methods). Sorting makes
% it easier to caluclate percentile p values
bootstrapped_perm_log_LASSO_left(:,1) = sort(mean(datasample(perm_left_ATL_log_LASSO,10000,2))');
bootstrapped_perm_log_LASSO_left(:,2) = sort(mean(datasample(perm_left_ATL_ant_log_LASSO,10000,2))');
bootstrapped_perm_log_LASSO_left(:,3) = sort(mean(datasample(perm_left_ATL_pos_log_LASSO,10000,2))');
bootstrapped_perm_log_LASSO_right(:,1) = sort(mean(datasample(perm_right_ATL_log_LASSO,10000,2))');
bootstrapped_perm_log_LASSO_right(:,2) = sort(mean(datasample(perm_right_ATL_ant_log_LASSO,10000,2))');
bootstrapped_perm_log_LASSO_right(:,3) = sort(mean(datasample(perm_right_ATL_pos_log_LASSO,10000,2))');
bootstrapped_perm_SOSLASSO_left(:,1) = sort(mean(datasample(perm_left_ATL_SOSLASSO,10000,2))');
bootstrapped_perm_SOSLASSO_left(:,2) = sort(mean(datasample(perm_left_ATL_ant_SOSLASSO,10000,2))');
bootstrapped_perm_SOSLASSO_left(:,3) = sort(mean(datasample(perm_left_ATL_pos_SOSLASSO,10000,2))');
bootstrapped_perm_SOSLASSO_right(:,1) = sort(mean(datasample(perm_right_ATL_SOSLASSO,10000,2))');
bootstrapped_perm_SOSLASSO_right(:,2) = sort(mean(datasample(perm_right_ATL_ant_SOSLASSO,10000,2))');
bootstrapped_perm_SOSLASSO_right(:,3) = sort(mean(datasample(perm_right_ATL_pos_SOSLASSO,10000,2))');

% do stats
classification_p_val = zeros(4,3);
m = mean(group_log_LASSO_left(:,1));
b = sum(bootstrapped_perm_log_LASSO_left(:,1) > m);
classification_p_val(1,1) = (b + 1)/(10000 + 1);
m = mean(group_log_LASSO_left(:,2));
b = sum(bootstrapped_perm_log_LASSO_left(:,2) > m);
classification_p_val(1,2) = (b + 1)/(10000 + 1);
m = mean(group_log_LASSO_left(:,3));
b = sum(bootstrapped_perm_log_LASSO_left(:,3) > m);
classification_p_val(1,3) = (b + 1)/(10000 + 1);
m = mean(group_log_LASSO_right(:,1));
b = sum(bootstrapped_perm_log_LASSO_right(:,1) > m);
classification_p_val(2,1) = (b + 1)/(10000 + 1);
m = mean(group_log_LASSO_right(:,2));
b = sum(bootstrapped_perm_log_LASSO_right(:,2) > m);
classification_p_val(2,2) = (b + 1)/(10000 + 1);
m = mean(group_log_LASSO_right(:,3));
b = sum(bootstrapped_perm_log_LASSO_right(:,3) > m);
classification_p_val(2,3) = (b + 1)/(10000 + 1);
m = mean(group_SOSLASSO_left(:,1));
b = sum(bootstrapped_perm_SOSLASSO_left(:,1) > m);
classification_p_val(3,1) = (b + 1)/(10000 + 1);
m = mean(group_SOSLASSO_left(:,2));
b = sum(bootstrapped_perm_SOSLASSO_left(:,2) > m);
classification_p_val(3,2) = (b + 1)/(10000 + 1);
m = mean(group_SOSLASSO_left(:,3));
b = sum(bootstrapped_perm_SOSLASSO_left(:,3) > m);
classification_p_val(3,3) = (b + 1)/(10000 + 1);
m = mean(group_SOSLASSO_right(:,1));
b = sum(bootstrapped_perm_SOSLASSO_right(:,1) > m);
classification_p_val(4,1) = (b + 1)/(10000 + 1);
m = mean(group_SOSLASSO_right(:,2));
b = sum(bootstrapped_perm_SOSLASSO_right(:,2) > m);
classification_p_val(4,2) = (b + 1)/(10000 + 1);
m = mean(group_SOSLASSO_right(:,3));
b = sum(bootstrapped_perm_SOSLASSO_right(:,3) > m);
classification_p_val(4,3) = (b + 1)/(10000 + 1);

% control the false discovery rate
for i = 1:size(classification_p_val,1)
    p = classification_p_val(i,:);
    classification_p_val(i,:) = mafdr(p,'BHFDR',1);
end

% test for differences between methods
classification_difference_p_val = zeros(4,3);
tmp = sort(mean(datasample((perm_left_ATL_log_LASSO - perm_left_ATL_SOSLASSO),10000,2))');
m = mean(group_log_LASSO_left(:,1) - group_SOSLASSO_left(:,1));
b = sum(tmp > m);
classification_difference_p_val(1,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_left_ATL_ant_log_LASSO - perm_left_ATL_ant_SOSLASSO),10000,2))');
m = mean(group_log_LASSO_left(:,2) - group_SOSLASSO_left(:,2));
b = sum(tmp > m);
classification_difference_p_val(1,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_left_ATL_pos_log_LASSO - perm_left_ATL_pos_SOSLASSO),10000,2))');
m = mean(group_log_LASSO_left(:,3) - group_SOSLASSO_left(:,3));
b = sum(tmp > m);
classification_difference_p_val(1,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_right_ATL_log_LASSO - perm_right_ATL_SOSLASSO),10000,2))');
m = mean(group_log_LASSO_right(:,1) - group_SOSLASSO_right(:,1));
b = sum(tmp > m);
classification_difference_p_val(2,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_right_ATL_ant_log_LASSO - perm_right_ATL_ant_SOSLASSO),10000,2))');
m = mean(group_log_LASSO_right(:,2) - group_SOSLASSO_right(:,2));
b = sum(tmp > m);
classification_difference_p_val(2,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_right_ATL_pos_log_LASSO - perm_right_ATL_pos_SOSLASSO),10000,2))');
m = mean(group_log_LASSO_right(:,3) - group_SOSLASSO_right(:,3));
b = sum(tmp > m);
classification_difference_p_val(2,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_left_ATL_SOSLASSO - perm_left_ATL_log_LASSO),10000,2))');
m = mean(group_SOSLASSO_left(:,1) - group_log_LASSO_left(:,1));
b = sum(tmp > m);
classification_difference_p_val(3,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_left_ATL_ant_SOSLASSO - perm_left_ATL_ant_log_LASSO),10000,2))');
m = mean(group_SOSLASSO_left(:,2) - group_log_LASSO_left(:,2));
b = sum(tmp > m);
classification_difference_p_val(3,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_left_ATL_ant_SOSLASSO - perm_left_ATL_ant_log_LASSO),10000,2))');
m = mean(group_SOSLASSO_left(:,2) - group_log_LASSO_left(:,2));
b = sum(tmp > m);
classification_difference_p_val(3,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_left_ATL_pos_SOSLASSO - perm_left_ATL_pos_log_LASSO),10000,2))');
m = mean(group_SOSLASSO_left(:,3) - group_log_LASSO_left(:,3));
b = sum(tmp > m);
classification_difference_p_val(3,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_right_ATL_SOSLASSO - perm_right_ATL_log_LASSO),10000,2))');
m = mean(group_SOSLASSO_right(:,1) - group_log_LASSO_right(:,1));
b = sum(tmp > m);
classification_difference_p_val(4,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_right_ATL_ant_SOSLASSO - perm_right_ATL_ant_log_LASSO),10000,2))');
m = mean(group_SOSLASSO_right(:,2) - group_log_LASSO_right(:,2));
b = sum(tmp > m);
classification_difference_p_val(4,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_right_ATL_pos_SOSLASSO - perm_right_ATL_pos_log_LASSO),10000,2))');
m = mean(group_SOSLASSO_right(:,3) - group_log_LASSO_right(:,3));
b = sum(tmp > m);
classification_difference_p_val(4,3) = (b + 1)/(10000 + 1);

% control the false discovery rate
for i = 1:size(classification_difference_p_val,1)
    p = classification_difference_p_val(i,:);
    classification_difference_p_val(i,:) = mafdr(p,'BHFDR',1);
end

% 1D. Plot results for the ATL

% plot
fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_log_LASSO_left));
b.FaceColor = 'flat';
b.CData = colours_LASSO_left(1:3:7,:);
set(gca,'xticklabel',{})
ylabel('Accuracy')
set(gca,'FontSize',20)
ylim([0.4,1])
yline(0.5,'--')
hold on
confint = 1.96*(std(group_log_LASSO_left)/sqrt(size(group_log_LASSO_left,1)));
errorbar(1:size(group_log_LASSO_left,2),mean(group_log_LASSO_left,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(classification_p_val,2)
    if (classification_p_val(1,i) < 0.05 && classification_difference_p_val(1,i) < 0.05)
       stars{i} = '*†';
    elseif classification_p_val(1,i) < 0.05
        stars{i} = '*';
    elseif classification_p_val(1,i) >= 0.05
       stars{i} = '';
    end
end
text((1:3)-0.1,mean(group_log_LASSO_left)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/log_LASSO_left.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_log_LASSO_right));
b.FaceColor = 'flat';
b.CData = colours_LASSO_right(1:3:7,:);
set(gca,'xticklabel',{})
ylabel('Accuracy')
set(gca,'FontSize',20)
ylim([0.4,1])
yline(0.5,'--')
hold on
confint = 1.96*(std(group_log_LASSO_right)/sqrt(size(group_log_LASSO_right,1)));
errorbar(1:size(group_log_LASSO_right,2),mean(group_log_LASSO_right,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(classification_p_val,2)
    if (classification_p_val(2,i) < 0.05 && classification_difference_p_val(2,i) < 0.05)
       stars{i} = '*†';
    elseif classification_p_val(2,i) < 0.05
        stars{i} = '*';
    elseif classification_p_val(2,i) >= 0.05
       stars{i} = '';
    end
end
% hard_coded fix for stars going out of the pane
tmp = [0.975,0.6707,0.975];
text((1:3)-0.1,tmp,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/log_LASSO_right.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_SOSLASSO_left));
b.FaceColor = 'flat';
b.CData = colours_SOSLASSO_left(1:3:7,:);
set(gca,'xticklabel',{})
ylabel('Accuracy')
set(gca,'FontSize',20)
ylim([0.4,1])
yline(0.5,'--')
hold on
confint = 1.96*(std(group_SOSLASSO_left)/sqrt(size(group_SOSLASSO_left,1)));
errorbar(1:size(group_SOSLASSO_left,2),mean(group_SOSLASSO_left,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(classification_p_val,2)
    if (classification_p_val(3,i) < 0.05 && classification_difference_p_val(3,i) < 0.05)
       stars{i} = '*†';
    elseif classification_p_val(3,i) < 0.05
        stars{i} = '*';
    elseif classification_p_val(3,i) >= 0.05
       stars{i} = '';
    end
end
% hard-coded fix for stars going out of the pane
tmp = [0.975,0.8467,0.975];
text((1:3)-0.1,tmp,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/SOSLASSO_left.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_SOSLASSO_right));
b.FaceColor = 'flat';
b.CData = colours_SOSLASSO_right(1:3:7,:);
set(gca,'xticklabel',{})
ylabel('Accuracy')
set(gca,'FontSize',20)
ylim([0.4,1])
yline(0.5,'--')
hold on
confint = 1.96*(std(group_SOSLASSO_right)/sqrt(size(group_SOSLASSO_right,1)));
errorbar(1:size(group_SOSLASSO_right,2),mean(group_SOSLASSO_right,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(classification_p_val,2)
    if (classification_p_val(4,i) < 0.05 && classification_difference_p_val(4,i) < 0.05)
       stars{i} = '*†';
    elseif classification_p_val(4,i) < 0.05
        stars{i} = '*';
    elseif classification_p_val(4,i) >= 0.05
       stars{i} = '';
    end
end
% hard_coded fix for stars going out of the pane
tmp = [0.975,0.6655,0.975];
text((1:3)-0.1,tmp,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/SOSLASSO_right.png')
close(fig)

% 2. CORRELATION (RSL with linear LASSO and with grOWL)

% get target dimensions
load('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/dilkina_norms.mat');
[C,z] = embed_similarity_matrix(dilkina_norms,3);
U = rescale_embedding(C,z);

% 2A. Organise linear LASSO results

whole_brain_lin_LASSO_D1 = zeros(10,size(subcode,2));
whole_brain_lin_LASSO_D2 = zeros(10,size(subcode,2));
whole_brain_lin_LASSO_D3 = zeros(10,size(subcode,2));
left_ATL_lin_LASSO_D1 = zeros(10,size(subcode,2));
left_ATL_lin_LASSO_D2 = zeros(10,size(subcode,2));
left_ATL_lin_LASSO_D3 = zeros(10,size(subcode,2));
right_ATL_lin_LASSO_D1 = zeros(10,size(subcode,2));
right_ATL_lin_LASSO_D2 = zeros(10,size(subcode,2));
right_ATL_lin_LASSO_D3 = zeros(10,size(subcode,2));
left_ATL_ant_lin_LASSO_D1 = zeros(10,size(subcode,2));
left_ATL_ant_lin_LASSO_D2 = zeros(10,size(subcode,2));
left_ATL_ant_lin_LASSO_D3 = zeros(10,size(subcode,2));
left_ATL_pos_lin_LASSO_D1 = zeros(10,size(subcode,2));
left_ATL_pos_lin_LASSO_D2 = zeros(10,size(subcode,2));
left_ATL_pos_lin_LASSO_D3 = zeros(10,size(subcode,2));
right_ATL_ant_lin_LASSO_D1 = zeros(10,size(subcode,2));
right_ATL_ant_lin_LASSO_D2 = zeros(10,size(subcode,2));
right_ATL_ant_lin_LASSO_D3 = zeros(10,size(subcode,2));
right_ATL_pos_lin_LASSO_D1 = zeros(10,size(subcode,2));
right_ATL_pos_lin_LASSO_D2 = zeros(10,size(subcode,2));
right_ATL_pos_lin_LASSO_D3 = zeros(10,size(subcode,2));
% perm results are arranged !! participant x random seed
perm_whole_brain_lin_LASSO_D1 = zeros(size(subcode,2),100);
perm_whole_brain_lin_LASSO_D2 = zeros(size(subcode,2),100);
perm_whole_brain_lin_LASSO_D3 = zeros(size(subcode,2),100);
perm_left_ATL_lin_LASSO_D1 = zeros(size(subcode,2),100);
perm_left_ATL_lin_LASSO_D2 = zeros(size(subcode,2),100);
perm_left_ATL_lin_LASSO_D3 = zeros(size(subcode,2),100);
perm_right_ATL_lin_LASSO_D1 = zeros(size(subcode,2),100);
perm_right_ATL_lin_LASSO_D2 = zeros(size(subcode,2),100);
perm_right_ATL_lin_LASSO_D3 = zeros(size(subcode,2),100);
perm_left_ATL_ant_lin_LASSO_D1 = zeros(size(subcode,2),100);
perm_left_ATL_ant_lin_LASSO_D2 = zeros(size(subcode,2),100);
perm_left_ATL_ant_lin_LASSO_D3 = zeros(size(subcode,2),100);
perm_left_ATL_pos_lin_LASSO_D1 = zeros(size(subcode,2),100);
perm_left_ATL_pos_lin_LASSO_D2 = zeros(size(subcode,2),100);
perm_left_ATL_pos_lin_LASSO_D3 = zeros(size(subcode,2),100);
perm_right_ATL_ant_lin_LASSO_D1 = zeros(size(subcode,2),100);
perm_right_ATL_ant_lin_LASSO_D2 = zeros(size(subcode,2),100);
perm_right_ATL_ant_lin_LASSO_D3 = zeros(size(subcode,2),100);
perm_right_ATL_pos_lin_LASSO_D1 = zeros(size(subcode,2),100);
perm_right_ATL_pos_lin_LASSO_D2 = zeros(size(subcode,2),100);
perm_right_ATL_pos_lin_LASSO_D3 = zeros(size(subcode,2),100);

whole_brain_lin_LASSO_D1_animate = zeros(10,size(subcode,2));
whole_brain_lin_LASSO_D2_animate = zeros(10,size(subcode,2));
whole_brain_lin_LASSO_D3_animate = zeros(10,size(subcode,2));
left_ATL_lin_LASSO_D1_animate = zeros(10,size(subcode,2));
left_ATL_lin_LASSO_D2_animate = zeros(10,size(subcode,2));
left_ATL_lin_LASSO_D3_animate = zeros(10,size(subcode,2));
right_ATL_lin_LASSO_D1_animate = zeros(10,size(subcode,2));
right_ATL_lin_LASSO_D2_animate = zeros(10,size(subcode,2));
right_ATL_lin_LASSO_D3_animate = zeros(10,size(subcode,2));
left_ATL_ant_lin_LASSO_D1_animate = zeros(10,size(subcode,2));
left_ATL_ant_lin_LASSO_D2_animate = zeros(10,size(subcode,2));
left_ATL_ant_lin_LASSO_D3_animate = zeros(10,size(subcode,2));
left_ATL_pos_lin_LASSO_D1_animate = zeros(10,size(subcode,2));
left_ATL_pos_lin_LASSO_D2_animate = zeros(10,size(subcode,2));
left_ATL_pos_lin_LASSO_D3_animate = zeros(10,size(subcode,2));
right_ATL_ant_lin_LASSO_D1_animate = zeros(10,size(subcode,2));
right_ATL_ant_lin_LASSO_D2_animate = zeros(10,size(subcode,2));
right_ATL_ant_lin_LASSO_D3_animate = zeros(10,size(subcode,2));
right_ATL_pos_lin_LASSO_D1_animate = zeros(10,size(subcode,2));
right_ATL_pos_lin_LASSO_D2_animate = zeros(10,size(subcode,2));
right_ATL_pos_lin_LASSO_D3_animate = zeros(10,size(subcode,2));
% perm results are arranged !! participant x random seed x fold (which will
% be averaged over
perm_whole_brain_lin_LASSO_D1_animate = zeros(size(subcode,2),100, 10);
perm_whole_brain_lin_LASSO_D2_animate = zeros(size(subcode,2),100, 10);
perm_whole_brain_lin_LASSO_D3_animate = zeros(size(subcode,2),100, 10);
perm_left_ATL_lin_LASSO_D1_animate = zeros(size(subcode,2),100, 10);
perm_left_ATL_lin_LASSO_D2_animate = zeros(size(subcode,2),100, 10);
perm_left_ATL_lin_LASSO_D3_animate = zeros(size(subcode,2),100, 10);
perm_right_ATL_lin_LASSO_D1_animate = zeros(size(subcode,2),100, 10);
perm_right_ATL_lin_LASSO_D2_animate = zeros(size(subcode,2),100, 10);
perm_right_ATL_lin_LASSO_D3_animate = zeros(size(subcode,2),100, 10);
perm_left_ATL_ant_lin_LASSO_D1_animate = zeros(size(subcode,2),100, 10);
perm_left_ATL_ant_lin_LASSO_D2_animate = zeros(size(subcode,2),100, 10);
perm_left_ATL_ant_lin_LASSO_D3_animate = zeros(size(subcode,2),100, 10);
perm_left_ATL_pos_lin_LASSO_D1_animate = zeros(size(subcode,2),100, 10);
perm_left_ATL_pos_lin_LASSO_D2_animate = zeros(size(subcode,2),100, 10);
perm_left_ATL_pos_lin_LASSO_D3_animate = zeros(size(subcode,2),100, 10);
perm_right_ATL_ant_lin_LASSO_D1_animate = zeros(size(subcode,2),100, 10);
perm_right_ATL_ant_lin_LASSO_D2_animate = zeros(size(subcode,2),100, 10);
perm_right_ATL_ant_lin_LASSO_D3_animate = zeros(size(subcode,2),100, 10);
perm_right_ATL_pos_lin_LASSO_D1_animate = zeros(size(subcode,2),100, 10);
perm_right_ATL_pos_lin_LASSO_D2_animate = zeros(size(subcode,2),100, 10);
perm_right_ATL_pos_lin_LASSO_D3_animate = zeros(size(subcode,2),100, 10);

whole_brain_lin_LASSO_D1_inanimate = zeros(10,size(subcode,2));
whole_brain_lin_LASSO_D2_inanimate = zeros(10,size(subcode,2));
whole_brain_lin_LASSO_D3_inanimate = zeros(10,size(subcode,2));
left_ATL_lin_LASSO_D1_inanimate = zeros(10,size(subcode,2));
left_ATL_lin_LASSO_D2_inanimate = zeros(10,size(subcode,2));
left_ATL_lin_LASSO_D3_inanimate = zeros(10,size(subcode,2));
right_ATL_lin_LASSO_D1_inanimate = zeros(10,size(subcode,2));
right_ATL_lin_LASSO_D2_inanimate = zeros(10,size(subcode,2));
right_ATL_lin_LASSO_D3_inanimate = zeros(10,size(subcode,2));
left_ATL_ant_lin_LASSO_D1_inanimate = zeros(10,size(subcode,2));
left_ATL_ant_lin_LASSO_D2_inanimate = zeros(10,size(subcode,2));
left_ATL_ant_lin_LASSO_D3_inanimate = zeros(10,size(subcode,2));
left_ATL_pos_lin_LASSO_D1_inanimate = zeros(10,size(subcode,2));
left_ATL_pos_lin_LASSO_D2_inanimate = zeros(10,size(subcode,2));
left_ATL_pos_lin_LASSO_D3_inanimate = zeros(10,size(subcode,2));
right_ATL_ant_lin_LASSO_D1_inanimate = zeros(10,size(subcode,2));
right_ATL_ant_lin_LASSO_D2_inanimate = zeros(10,size(subcode,2));
right_ATL_ant_lin_LASSO_D3_inanimate = zeros(10,size(subcode,2));
right_ATL_pos_lin_LASSO_D1_inanimate = zeros(10,size(subcode,2));
right_ATL_pos_lin_LASSO_D2_inanimate = zeros(10,size(subcode,2));
right_ATL_pos_lin_LASSO_D3_inanimate = zeros(10,size(subcode,2));
perm_whole_brain_lin_LASSO_D1_inanimate = zeros(size(subcode,2),100, 10);
perm_whole_brain_lin_LASSO_D2_inanimate = zeros(size(subcode,2),100, 10);
perm_whole_brain_lin_LASSO_D3_inanimate = zeros(size(subcode,2),100, 10);
perm_left_ATL_lin_LASSO_D1_inanimate = zeros(size(subcode,2),100, 10);
perm_left_ATL_lin_LASSO_D2_inanimate = zeros(size(subcode,2),100, 10);
perm_left_ATL_lin_LASSO_D3_inanimate = zeros(size(subcode,2),100, 10);
perm_right_ATL_lin_LASSO_D1_inanimate = zeros(size(subcode,2),100, 10);
perm_right_ATL_lin_LASSO_D2_inanimate = zeros(size(subcode,2),100, 10);
perm_right_ATL_lin_LASSO_D3_inanimate = zeros(size(subcode,2),100, 10);
perm_left_ATL_ant_lin_LASSO_D1_inanimate = zeros(size(subcode,2),100, 10);
perm_left_ATL_ant_lin_LASSO_D2_inanimate = zeros(size(subcode,2),100, 10);
perm_left_ATL_ant_lin_LASSO_D3_inanimate = zeros(size(subcode,2),100, 10);
perm_left_ATL_pos_lin_LASSO_D1_inanimate = zeros(size(subcode,2),100, 10);
perm_left_ATL_pos_lin_LASSO_D2_inanimate = zeros(size(subcode,2),100, 10);
perm_left_ATL_pos_lin_LASSO_D3_inanimate = zeros(size(subcode,2),100, 10);
perm_right_ATL_ant_lin_LASSO_D1_inanimate = zeros(size(subcode,2),100, 10);
perm_right_ATL_ant_lin_LASSO_D2_inanimate = zeros(size(subcode,2),100, 10);
perm_right_ATL_ant_lin_LASSO_D3_inanimate = zeros(size(subcode,2),100, 10);
perm_right_ATL_pos_lin_LASSO_D1_inanimate = zeros(size(subcode,2),100, 10);
perm_right_ATL_pos_lin_LASSO_D2_inanimate = zeros(size(subcode,2),100, 10);
perm_right_ATL_pos_lin_LASSO_D3_inanimate = zeros(size(subcode,2),100, 10);

% for every participant
for s = 1:size(subcode,2)

    % load results
    result = load([root,'/LASSO/linear/tedana/',subcode{s},'/results.mat']);

    % collate 
    whole_brain_lin_LASSO_D1(:,s) = result.averaged_whole_brain_D1.output;
    whole_brain_lin_LASSO_D2(:,s) = result.averaged_whole_brain_D2.output;
    whole_brain_lin_LASSO_D3(:,s) = result.averaged_whole_brain_D3.output;
    left_ATL_lin_LASSO_D1(:,s) = result.averaged_left_ATL_D1.output;
    left_ATL_lin_LASSO_D2(:,s) = result.averaged_left_ATL_D2.output;
    left_ATL_lin_LASSO_D3(:,s) = result.averaged_left_ATL_D3.output;
    right_ATL_lin_LASSO_D1(:,s) = result.averaged_right_ATL_D1.output;
    right_ATL_lin_LASSO_D2(:,s) = result.averaged_right_ATL_D2.output;
    right_ATL_lin_LASSO_D3(:,s) = result.averaged_right_ATL_D3.output;
    left_ATL_ant_lin_LASSO_D1(:,s) = result.averaged_left_ATL_ant_D1.output;
    left_ATL_ant_lin_LASSO_D2(:,s) = result.averaged_left_ATL_ant_D2.output;
    left_ATL_ant_lin_LASSO_D3(:,s) = result.averaged_left_ATL_ant_D3.output;
    left_ATL_pos_lin_LASSO_D1(:,s) = result.averaged_left_ATL_pos_D1.output;
    left_ATL_pos_lin_LASSO_D2(:,s) = result.averaged_left_ATL_pos_D2.output;
    left_ATL_pos_lin_LASSO_D3(:,s) = result.averaged_left_ATL_pos_D3.output;
    right_ATL_ant_lin_LASSO_D1(:,s) = result.averaged_right_ATL_ant_D1.output;
    right_ATL_ant_lin_LASSO_D2(:,s) = result.averaged_right_ATL_ant_D2.output;
    right_ATL_ant_lin_LASSO_D3(:,s) = result.averaged_right_ATL_ant_D3.output;
    right_ATL_pos_lin_LASSO_D1(:,s) = result.averaged_right_ATL_pos_D1.output;
    right_ATL_pos_lin_LASSO_D2(:,s) = result.averaged_right_ATL_pos_D2.output;
    right_ATL_pos_lin_LASSO_D3(:,s) = result.averaged_right_ATL_pos_D3.output;

    perm_whole_brain_lin_LASSO_D1(s,:) = mean(result.averaged_whole_brain_D1.permoutput)';
    perm_whole_brain_lin_LASSO_D2(s,:) = mean(result.averaged_whole_brain_D2.permoutput)';
    perm_whole_brain_lin_LASSO_D3(s,:) = mean(result.averaged_whole_brain_D3.permoutput)';
    perm_left_ATL_lin_LASSO_D1(s,:) = mean(result.averaged_left_ATL_D1.permoutput)';
    perm_left_ATL_lin_LASSO_D2(s,:) = mean(result.averaged_left_ATL_D2.permoutput)';
    perm_left_ATL_lin_LASSO_D3(s,:) = mean(result.averaged_left_ATL_D3.permoutput)';
    perm_right_ATL_lin_LASSO_D1(s,:) = mean(result.averaged_right_ATL_D1.permoutput)';
    perm_right_ATL_lin_LASSO_D2(s,:) = mean(result.averaged_right_ATL_D2.permoutput)';
    perm_right_ATL_lin_LASSO_D3(s,:) = mean(result.averaged_right_ATL_D3.permoutput)';
    perm_left_ATL_ant_lin_LASSO_D1(s,:) = mean(result.averaged_left_ATL_ant_D1.permoutput)';
    perm_left_ATL_ant_lin_LASSO_D2(s,:) = mean(result.averaged_left_ATL_ant_D2.permoutput)';
    perm_left_ATL_ant_lin_LASSO_D3(s,:) = mean(result.averaged_left_ATL_ant_D3.permoutput)';
    perm_left_ATL_pos_lin_LASSO_D1(s,:) = mean(result.averaged_left_ATL_pos_D1.permoutput)';
    perm_left_ATL_pos_lin_LASSO_D2(s,:) = mean(result.averaged_left_ATL_pos_D2.permoutput)';
    perm_left_ATL_pos_lin_LASSO_D3(s,:) = mean(result.averaged_left_ATL_pos_D3.permoutput)';
    perm_right_ATL_ant_lin_LASSO_D1(s,:) = mean(result.averaged_right_ATL_ant_D1.permoutput)';
    perm_right_ATL_ant_lin_LASSO_D2(s,:) = mean(result.averaged_right_ATL_ant_D2.permoutput)';
    perm_right_ATL_ant_lin_LASSO_D3(s,:) = mean(result.averaged_right_ATL_ant_D3.permoutput)';
    perm_right_ATL_pos_lin_LASSO_D1(s,:) = mean(result.averaged_right_ATL_pos_D1.permoutput)';
    perm_right_ATL_pos_lin_LASSO_D2(s,:) = mean(result.averaged_right_ATL_pos_D2.permoutput)';
    perm_right_ATL_pos_lin_LASSO_D3(s,:) = mean(result.averaged_right_ATL_pos_D3.permoutput)';

    % calculate within-domain correlations
    % for each holdout fold
    for ho = 1:10 
        % dimension 1
            % get target for that holdout set
            utst = U(metadata(s).cvind(:,1)==ho,1);
            % WHOLE BRAIN
            % get predicted coordinates
            allpred = result.averaged_whole_brain_D1.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            whole_brain_lin_LASSO_D1_animate(ho,s) = corr(pred(1:5),utst(1:5));
            whole_brain_lin_LASSO_D1_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                whole_brain_lin_LASSO_D1_animate(ho,s) = 0;
                whole_brain_lin_LASSO_D1_inanimate(ho,s) = 0;
            end
            % LEFT ATL
            % get predicted coordinates
            allpred = result.averaged_left_ATL_D1.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            left_ATL_lin_LASSO_D1_animate(ho,s) = corr(pred(1:5),utst(1:5));
            left_ATL_lin_LASSO_D1_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                left_ATL_lin_LASSO_D1_animate(ho,s) = 0;
                left_ATL_lin_LASSO_D1_inanimate(ho,s) = 0;
            end
            % RIGHT ATL
            allpred = result.averaged_right_ATL_D1.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            right_ATL_lin_LASSO_D1_animate(ho,s) = corr(pred(1:5),utst(1:5));
            right_ATL_lin_LASSO_D1_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                right_ATL_lin_LASSO_D1_animate(ho,s) = 0;
                right_ATL_lin_LASSO_D1_inanimate(ho,s) = 0;
            end
            % LEFT ATL (ANTERIOR)
            % get predicted coordinates
            allpred = result.averaged_left_ATL_ant_D1.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            left_ATL_ant_lin_LASSO_D1_animate(ho,s) = corr(pred(1:5),utst(1:5));
            left_ATL_ant_lin_LASSO_D1_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                left_ATL_ant_lin_LASSO_D1_animate(ho,s) = 0;
                left_ATL_ant_lin_LASSO_D1_inanimate(ho,s) = 0;
            end
            % LEFT ATL (POSTERIOR)
            % get predicted coordinates
            allpred = result.averaged_left_ATL_pos_D1.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            left_ATL_pos_lin_LASSO_D1_animate(ho,s) = corr(pred(1:5),utst(1:5));
            left_ATL_pos_lin_LASSO_D1_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                left_ATL_pos_lin_LASSO_D1_animate(ho,s) = 0;
                left_ATL_pos_lin_LASSO_D1_inanimate(ho,s) = 0;
            end
            % RIGHT ATL (ANTERIOR)
            % get predicted coordinates
            allpred = result.averaged_right_ATL_ant_D1.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            right_ATL_ant_lin_LASSO_D1_animate(ho,s) = corr(pred(1:5),utst(1:5));
            right_ATL_ant_lin_LASSO_D1_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                right_ATL_ant_lin_LASSO_D1_animate(ho,s) = 0;
                right_ATL_ant_lin_LASSO_D1_inanimate(ho,s) = 0;
            end
            % RIGHT ATL (POSTERIOR)
            % get predicted coordinates
            allpred = result.averaged_right_ATL_pos_D1.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            right_ATL_pos_lin_LASSO_D1_animate(ho,s) = corr(pred(1:5),utst(1:5));
            right_ATL_pos_lin_LASSO_D1_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                right_ATL_pos_lin_LASSO_D1_animate(ho,s) = 0;
                right_ATL_pos_lin_LASSO_D1_inanimate(ho,s) = 0;
            end

        % Dimension 2
            % get target for that holdout set
            utst = U(metadata(s).cvind(:,1)==ho,2);
            % WHOLE BRAIN
            % get predicted coordinates
            allpred = result.averaged_whole_brain_D2.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            whole_brain_lin_LASSO_D2_animate(ho,s) = corr(pred(1:5),utst(1:5));
            whole_brain_lin_LASSO_D2_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                whole_brain_lin_LASSO_D2_animate(ho,s) = 0;
                whole_brain_lin_LASSO_D2_inanimate(ho,s) = 0;
            end
            % LEFT ATL
            % get predicted coordinates
            allpred = result.averaged_left_ATL_D2.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            left_ATL_lin_LASSO_D2_animate(ho,s) = corr(pred(1:5),utst(1:5));
            left_ATL_lin_LASSO_D2_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                left_ATL_lin_LASSO_D2_animate(ho,s) = 0;
                left_ATL_lin_LASSO_D2_inanimate(ho,s) = 0;
            end
            % RIGHT ATL
            allpred = result.averaged_right_ATL_D2.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            right_ATL_lin_LASSO_D2_animate(ho,s) = corr(pred(1:5),utst(1:5));
            right_ATL_lin_LASSO_D2_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                right_ATL_lin_LASSO_D2_animate(ho,s) = 0;
                right_ATL_lin_LASSO_D2_inanimate(ho,s) = 0;
            end
            % LEFT ATL (ANTERIOR)
            % get predicted coordinates
            allpred = result.averaged_left_ATL_ant_D2.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            left_ATL_ant_lin_LASSO_D2_animate(ho,s) = corr(pred(1:5),utst(1:5));
            left_ATL_ant_lin_LASSO_D2_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                left_ATL_ant_lin_LASSO_D2_animate(ho,s) = 0;
                left_ATL_ant_lin_LASSO_D2_inanimate(ho,s) = 0;
            end
            % LEFT ATL (POSTERIOR)
            % get predicted coordinates
            allpred = result.averaged_left_ATL_pos_D2.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            left_ATL_pos_lin_LASSO_D2_animate(ho,s) = corr(pred(1:5),utst(1:5));
            left_ATL_pos_lin_LASSO_D2_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                left_ATL_pos_lin_LASSO_D2_animate(ho,s) = 0;
                left_ATL_pos_lin_LASSO_D2_inanimate(ho,s) = 0;
            end
            % RIGHT ATL (ANTERIOR)
            % get predicted coordinates
            allpred = result.averaged_right_ATL_ant_D2.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            right_ATL_ant_lin_LASSO_D2_animate(ho,s) = corr(pred(1:5),utst(1:5));
            right_ATL_ant_lin_LASSO_D2_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                right_ATL_ant_lin_LASSO_D2_animate(ho,s) = 0;
                right_ATL_ant_lin_LASSO_D2_inanimate(ho,s) = 0;
            end
            % RIGHT ATL (POSTERIOR)
            % get predicted coordinates
            allpred = result.averaged_right_ATL_pos_D2.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            right_ATL_pos_lin_LASSO_D2_animate(ho,s) = corr(pred(1:5),utst(1:5));
            right_ATL_pos_lin_LASSO_D2_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                right_ATL_pos_lin_LASSO_D2_animate(ho,s) = 0;
                right_ATL_pos_lin_LASSO_D2_inanimate(ho,s) = 0;
            end

        % Dimension 3
            % get target for that holdout set
            utst = U(metadata(s).cvind(:,1)==ho,3);
            % WHOLE BRAIN
            % get predicted coordinates
            allpred = result.averaged_whole_brain_D3.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            whole_brain_lin_LASSO_D3_animate(ho,s) = corr(pred(1:5),utst(1:5));
            whole_brain_lin_LASSO_D3_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                whole_brain_lin_LASSO_D3_animate(ho,s) = 0;
                whole_brain_lin_LASSO_D3_inanimate(ho,s) = 0;
            end
            % LEFT ATL
            % get predicted coordinates
            allpred = result.averaged_left_ATL_D3.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            left_ATL_lin_LASSO_D3_animate(ho,s) = corr(pred(1:5),utst(1:5));
            left_ATL_lin_LASSO_D3_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                left_ATL_lin_LASSO_D3_animate(ho,s) = 0;
                left_ATL_lin_LASSO_D3_inanimate(ho,s) = 0;
            end
            % RIGHT ATL
            allpred = result.averaged_right_ATL_D3.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            right_ATL_lin_LASSO_D3_animate(ho,s) = corr(pred(1:5),utst(1:5));
            right_ATL_lin_LASSO_D3_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                right_ATL_lin_LASSO_D3_animate(ho,s) = 0;
                right_ATL_lin_LASSO_D3_inanimate(ho,s) = 0;
            end
            % LEFT ATL (ANTERIOR)
            % get predicted coordinates
            allpred = result.averaged_left_ATL_ant_D3.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            left_ATL_ant_lin_LASSO_D3_animate(ho,s) = corr(pred(1:5),utst(1:5));
            left_ATL_ant_lin_LASSO_D3_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                left_ATL_ant_lin_LASSO_D3_animate(ho,s) = 0;
                left_ATL_ant_lin_LASSO_D3_inanimate(ho,s) = 0;
            end
            % LEFT ATL (POSTERIOR)
            % get predicted coordinates
            allpred = result.averaged_left_ATL_pos_D3.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            left_ATL_pos_lin_LASSO_D3_animate(ho,s) = corr(pred(1:5),utst(1:5));
            left_ATL_pos_lin_LASSO_D3_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                left_ATL_pos_lin_LASSO_D3_animate(ho,s) = 0;
                left_ATL_pos_lin_LASSO_D3_inanimate(ho,s) = 0;
            end
            % RIGHT ATL (ANTERIOR)
            % get predicted coordinates
            allpred = result.averaged_right_ATL_ant_D3.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            right_ATL_ant_lin_LASSO_D3_animate(ho,s) = corr(pred(1:5),utst(1:5));
            right_ATL_ant_lin_LASSO_D3_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                right_ATL_ant_lin_LASSO_D3_animate(ho,s) = 0;
                right_ATL_ant_lin_LASSO_D3_inanimate(ho,s) = 0;
            end
            % RIGHT ATL (POSTERIOR)
            % get predicted coordinates
            allpred = result.averaged_right_ATL_pos_D3.predictedcoords;
            % index just coordinates in that holdout set
            pred = allpred(metadata(s).cvind(:,1)==ho);
            % calculate within-domain correlation
            right_ATL_pos_lin_LASSO_D3_animate(ho,s) = corr(pred(1:5),utst(1:5));
            right_ATL_pos_lin_LASSO_D3_inanimate(ho,s)= corr(pred(6:10),utst(6:10));
            if var(pred) == 0
                right_ATL_pos_lin_LASSO_D3_animate(ho,s) = 0;
                right_ATL_pos_lin_LASSO_D3_inanimate(ho,s) = 0;
            end

            % Within-domain permutation correlations

            for seed = 1:100
                % dimension 1
                % get target for that holdout set
                utst = U(metadata(s).cvind(:,1)==ho,1);
                    % WHOLE BRAIN
                    % get predicted coordinates
                    allpred = result.averaged_whole_brain_D1.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_whole_brain_lin_LASSO_D1_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_whole_brain_lin_LASSO_D1_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_whole_brain_lin_LASSO_D1_animate(s,seed,ho) = 0;
                        perm_whole_brain_lin_LASSO_D1_inanimate(s,seed,ho) = 0;
                    end
                    % LEFT ATL
                    % get predicted coordinates
                    allpred = result.averaged_left_ATL_D1.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_left_ATL_lin_LASSO_D1_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_left_ATL_lin_LASSO_D1_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_left_ATL_lin_LASSO_D1_animate(s,seed,ho) = 0;
                        perm_left_ATL_lin_LASSO_D1_inanimate(s,seed,ho) = 0;
                    end
                    % RIGHT ATL
                    % get predicted coordinates
                    allpred = result.averaged_right_ATL_D1.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_right_ATL_lin_LASSO_D1_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_right_ATL_lin_LASSO_D1_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_right_ATL_lin_LASSO_D1_animate(s,seed,ho) = 0;
                        perm_right_ATL_lin_LASSO_D1_inanimate(s,seed,ho) = 0;
                    end
                    % LEFT ATL (ANTERIOR)
                    % get predicted coordinates
                    allpred = result.averaged_left_ATL_ant_D1.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_left_ATL_ant_lin_LASSO_D1_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_left_ATL_ant_lin_LASSO_D1_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_left_ATL_ant_lin_LASSO_D1_animate(s,seed,ho) = 0;
                        perm_left_ATL_ant_lin_LASSO_D1_inanimate(s,seed,ho) = 0;
                    end
                    % LEFT ATL (POSTERIOR)
                    % get predicted coordinates
                    allpred = result.averaged_left_ATL_pos_D1.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_left_ATL_pos_lin_LASSO_D1_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_left_ATL_pos_lin_LASSO_D1_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_left_ATL_pos_lin_LASSO_D1_animate(s,seed,ho) = 0;
                        perm_left_ATL_pos_lin_LASSO_D1_inanimate(s,seed,ho) = 0;
                    end
                    % RiGHT ATL (ANTERIOR)
                    % get predicted coordinates
                    allpred = result.averaged_right_ATL_ant_D1.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_right_ATL_ant_lin_LASSO_D1_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_right_ATL_ant_lin_LASSO_D1_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_right_ATL_ant_lin_LASSO_D1_animate(s,seed,ho) = 0;
                        perm_right_ATL_ant_lin_LASSO_D1_inanimate(s,seed,ho) = 0;
                    end
                    % RIGHT ATL (POSTERIOR)
                    % get predicted coordinates
                    allpred = result.averaged_right_ATL_pos_D1.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_right_ATL_pos_lin_LASSO_D1_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_right_ATL_pos_lin_LASSO_D1_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_right_ATL_pos_lin_LASSO_D1_animate(s,seed,ho) = 0;
                        perm_right_ATL_pos_lin_LASSO_D1_inanimate(s,seed,ho) = 0;
                    end

                % dimension 2
                % get target for that holdout set
                utst = U(metadata(s).cvind(:,1)==ho,2);
                    % WHOLE BRAIN
                    % get predicted coordinates
                    allpred = result.averaged_whole_brain_D2.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_whole_brain_lin_LASSO_D2_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_whole_brain_lin_LASSO_D2_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_whole_brain_lin_LASSO_D2_animate(s,seed,ho) = 0;
                        perm_whole_brain_lin_LASSO_D2_inanimate(s,seed,ho) = 0;
                    end
                    % LEFT ATL
                    % get predicted coordinates
                    allpred = result.averaged_left_ATL_D2.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_left_ATL_lin_LASSO_D2_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_left_ATL_lin_LASSO_D2_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_left_ATL_lin_LASSO_D2_animate(s,seed,ho) = 0;
                        perm_left_ATL_lin_LASSO_D2_inanimate(s,seed,ho) = 0;
                    end
                    % RIGHT ATL
                    % get predicted coordinates
                    allpred = result.averaged_right_ATL_D2.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_right_ATL_lin_LASSO_D2_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_right_ATL_lin_LASSO_D2_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_right_ATL_lin_LASSO_D2_animate(s,seed,ho) = 0;
                        perm_right_ATL_lin_LASSO_D2_inanimate(s,seed,ho) = 0;
                    end
                    % LEFT ATL (ANTERIOR)
                    % get predicted coordinates
                    allpred = result.averaged_left_ATL_ant_D2.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_left_ATL_ant_lin_LASSO_D2_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_left_ATL_ant_lin_LASSO_D2_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_left_ATL_ant_lin_LASSO_D2_animate(s,seed,ho) = 0;
                        perm_left_ATL_ant_lin_LASSO_D2_inanimate(s,seed,ho) = 0;
                    end
                    % LEFT ATL (POSTERIOR)
                    % get predicted coordinates
                    allpred = result.averaged_left_ATL_pos_D2.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_left_ATL_pos_lin_LASSO_D2_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_left_ATL_pos_lin_LASSO_D2_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_left_ATL_pos_lin_LASSO_D2_animate(s,seed,ho) = 0;
                        perm_left_ATL_pos_lin_LASSO_D2_inanimate(s,seed,ho) = 0;
                    end
                    % RiGHT ATL (ANTERIOR)
                    % get predicted coordinates
                    allpred = result.averaged_right_ATL_ant_D2.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_right_ATL_ant_lin_LASSO_D2_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_right_ATL_ant_lin_LASSO_D2_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_right_ATL_ant_lin_LASSO_D2_animate(s,seed,ho) = 0;
                        perm_right_ATL_ant_lin_LASSO_D2_inanimate(s,seed,ho) = 0;
                    end
                    % RIGHT ATL (POSTERIOR)
                    % get predicted coordinates
                    allpred = result.averaged_right_ATL_pos_D2.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_right_ATL_pos_lin_LASSO_D2_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_right_ATL_pos_lin_LASSO_D2_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_right_ATL_pos_lin_LASSO_D2_animate(s,seed,ho) = 0;
                        perm_right_ATL_pos_lin_LASSO_D2_inanimate(s,seed,ho) = 0;
                    end
                % dimension 3
                % get target for that holdout set
                utst = U(metadata(s).cvind(:,1)==ho,3);
                    % WHOLE BRAIN
                    % get predicted coordinates
                    allpred = result.averaged_whole_brain_D3.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_whole_brain_lin_LASSO_D3_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_whole_brain_lin_LASSO_D3_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_whole_brain_lin_LASSO_D3_animate(s,seed,ho) = 0;
                        perm_whole_brain_lin_LASSO_D3_inanimate(s,seed,ho) = 0;
                    end
                    % LEFT ATL
                    % get predicted coordinates
                    allpred = result.averaged_left_ATL_D3.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_left_ATL_lin_LASSO_D3_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_left_ATL_lin_LASSO_D3_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_left_ATL_lin_LASSO_D3_animate(s,seed,ho) = 0;
                        perm_left_ATL_lin_LASSO_D3_inanimate(s,seed,ho) = 0;
                    end
                    % RIGHT ATL
                    % get predicted coordinates
                    allpred = result.averaged_right_ATL_D3.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_right_ATL_lin_LASSO_D3_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_right_ATL_lin_LASSO_D3_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_right_ATL_lin_LASSO_D3_animate(s,seed,ho) = 0;
                        perm_right_ATL_lin_LASSO_D3_inanimate(s,seed,ho) = 0;
                    end
                    % LEFT ATL (ANTERIOR)
                    % get predicted coordinates
                    allpred = result.averaged_left_ATL_ant_D3.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_left_ATL_ant_lin_LASSO_D3_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_left_ATL_ant_lin_LASSO_D3_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_left_ATL_ant_lin_LASSO_D3_animate(s,seed,ho) = 0;
                        perm_left_ATL_ant_lin_LASSO_D3_inanimate(s,seed,ho) = 0;
                    end
                    % LEFT ATL (POSTERIOR)
                    % get predicted coordinates
                    allpred = result.averaged_left_ATL_pos_D3.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_left_ATL_pos_lin_LASSO_D3_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_left_ATL_pos_lin_LASSO_D3_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_left_ATL_pos_lin_LASSO_D3_animate(s,seed,ho) = 0;
                        perm_left_ATL_pos_lin_LASSO_D3_inanimate(s,seed,ho) = 0;
                    end
                    % RiGHT ATL (ANTERIOR)
                    % get predicted coordinates
                    allpred = result.averaged_right_ATL_ant_D3.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_right_ATL_ant_lin_LASSO_D3_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_right_ATL_ant_lin_LASSO_D3_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_right_ATL_ant_lin_LASSO_D3_animate(s,seed,ho) = 0;
                        perm_right_ATL_ant_lin_LASSO_D3_inanimate(s,seed,ho) = 0;
                    end
                    % RIGHT ATL (POSTERIOR)
                    % get predicted coordinates
                    allpred = result.averaged_right_ATL_pos_D3.permpredictedcoords(:,seed);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate within-domain correlation
                    perm_right_ATL_pos_lin_LASSO_D3_animate(s,seed,ho) = corr(pred(1:5),utst(1:5));
                    perm_right_ATL_pos_lin_LASSO_D3_inanimate(s,seed,ho)= corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        perm_right_ATL_pos_lin_LASSO_D3_animate(s,seed,ho) = 0;
                        perm_right_ATL_pos_lin_LASSO_D3_inanimate(s,seed,ho) = 0;
                    end
            end
    end
end

perm_whole_brain_lin_LASSO_D1_animate = mean(perm_whole_brain_lin_LASSO_D1_animate,3);
perm_whole_brain_lin_LASSO_D2_animate = mean(perm_whole_brain_lin_LASSO_D2_animate,3);
perm_whole_brain_lin_LASSO_D3_animate = mean(perm_whole_brain_lin_LASSO_D3_animate,3);
perm_left_ATL_lin_LASSO_D1_animate = mean(perm_left_ATL_lin_LASSO_D1_animate,3);
perm_left_ATL_lin_LASSO_D2_animate = mean(perm_left_ATL_lin_LASSO_D2_animate,3);
perm_left_ATL_lin_LASSO_D3_animate = mean(perm_left_ATL_lin_LASSO_D3_animate,3);
perm_right_ATL_lin_LASSO_D1_animate = mean(perm_right_ATL_lin_LASSO_D1_animate,3);
perm_right_ATL_lin_LASSO_D2_animate = mean(perm_right_ATL_lin_LASSO_D2_animate,3);
perm_right_ATL_lin_LASSO_D3_animate = mean(perm_right_ATL_lin_LASSO_D3_animate,3);
perm_left_ATL_ant_lin_LASSO_D1_animate = mean(perm_left_ATL_ant_lin_LASSO_D1_animate ,3);
perm_left_ATL_ant_lin_LASSO_D2_animate = mean(perm_left_ATL_ant_lin_LASSO_D2_animate,3);
perm_left_ATL_ant_lin_LASSO_D3_animate = mean(perm_left_ATL_ant_lin_LASSO_D3_animate,3);
perm_left_ATL_pos_lin_LASSO_D1_animate = mean(perm_left_ATL_pos_lin_LASSO_D1_animate,3);
perm_left_ATL_pos_lin_LASSO_D2_animate = mean(perm_left_ATL_pos_lin_LASSO_D2_animate,3);
perm_left_ATL_pos_lin_LASSO_D3_animate = mean(perm_left_ATL_pos_lin_LASSO_D3_animate,3);
perm_right_ATL_ant_lin_LASSO_D1_animate = mean(perm_right_ATL_ant_lin_LASSO_D1_animate,3);
perm_right_ATL_ant_lin_LASSO_D2_animate = mean(perm_right_ATL_ant_lin_LASSO_D2_animate,3);
perm_right_ATL_ant_lin_LASSO_D3_animate = mean(perm_right_ATL_ant_lin_LASSO_D3_animate,3);
perm_right_ATL_pos_lin_LASSO_D1_animate = mean(perm_right_ATL_pos_lin_LASSO_D1_animate,3);
perm_right_ATL_pos_lin_LASSO_D2_animate = mean(perm_right_ATL_pos_lin_LASSO_D2_animate,3);
perm_right_ATL_pos_lin_LASSO_D3_animate = mean(perm_right_ATL_pos_lin_LASSO_D3_animate,3);

perm_whole_brain_lin_LASSO_D1_inanimate = mean(perm_whole_brain_lin_LASSO_D1_inanimate,3);
perm_whole_brain_lin_LASSO_D2_inanimate = mean(perm_whole_brain_lin_LASSO_D2_inanimate,3);
perm_whole_brain_lin_LASSO_D3_inanimate = mean(perm_whole_brain_lin_LASSO_D3_inanimate,3);
perm_left_ATL_lin_LASSO_D1_inanimate = mean(perm_left_ATL_lin_LASSO_D1_inanimate,3);
perm_left_ATL_lin_LASSO_D2_inanimate = mean(perm_left_ATL_lin_LASSO_D2_inanimate,3);
perm_left_ATL_lin_LASSO_D3_inanimate = mean(perm_left_ATL_lin_LASSO_D3_inanimate,3);
perm_right_ATL_lin_LASSO_D1_inanimate = mean(perm_right_ATL_lin_LASSO_D1_inanimate,3);
perm_right_ATL_lin_LASSO_D2_inanimate = mean(perm_right_ATL_lin_LASSO_D2_inanimate,3);
perm_right_ATL_lin_LASSO_D3_inanimate = mean(perm_right_ATL_lin_LASSO_D3_inanimate,3);
perm_left_ATL_ant_lin_LASSO_D1_inanimate = mean(perm_left_ATL_ant_lin_LASSO_D1_inanimate ,3);
perm_left_ATL_ant_lin_LASSO_D2_inanimate = mean(perm_left_ATL_ant_lin_LASSO_D2_inanimate,3);
perm_left_ATL_ant_lin_LASSO_D3_inanimate = mean(perm_left_ATL_ant_lin_LASSO_D3_inanimate,3);
perm_left_ATL_pos_lin_LASSO_D1_inanimate = mean(perm_left_ATL_pos_lin_LASSO_D1_inanimate,3);
perm_left_ATL_pos_lin_LASSO_D2_inanimate = mean(perm_left_ATL_pos_lin_LASSO_D2_inanimate,3);
perm_left_ATL_pos_lin_LASSO_D3_inanimate = mean(perm_left_ATL_pos_lin_LASSO_D3_inanimate,3);
perm_right_ATL_ant_lin_LASSO_D1_inanimate = mean(perm_right_ATL_ant_lin_LASSO_D1_inanimate,3);
perm_right_ATL_ant_lin_LASSO_D2_inanimate = mean(perm_right_ATL_ant_lin_LASSO_D2_inanimate,3);
perm_right_ATL_ant_lin_LASSO_D3_inanimate = mean(perm_right_ATL_ant_lin_LASSO_D3_inanimate,3);
perm_right_ATL_pos_lin_LASSO_D1_inanimate = mean(perm_right_ATL_pos_lin_LASSO_D1_inanimate,3);
perm_right_ATL_pos_lin_LASSO_D2_inanimate = mean(perm_right_ATL_pos_lin_LASSO_D2_inanimate,3);
perm_right_ATL_pos_lin_LASSO_D3_inanimate = mean(perm_right_ATL_pos_lin_LASSO_D3_inanimate,3);

% 2B. Organise grOWL results

% initialise temporary results
result = {};
permresult = {};

% load results from condor

% for each ROI
for r = 1:length(ROI)
    
    % load results
    load([root,'/condor/derivatives/correlation/grOWL/performance/',ROI{r},'/final_performance.mat']);
    
    % initialise output
    output = zeros(10,size(subcode,2),3);
    % initialise within-domain output
    animateoutput = zeros(10,size(subcode,2),3);
    inanimateoutput = zeros(10,size(subcode,2),3);

    % for each participant
    for s = 1:size(subcode,2)
        % get just the rows applying to that participant
        tmp = str2num(erase(subcode{s},'sub-'));
        subtable = Tallcv(Tallcv.subject==tmp,:);
        
        % for each holdout fold
        for ho = 1:10 
            % for each dimension
            for d = 1:3
                % get target for that holdout set
                utst = U(metadata(s).cvind(:,1)==ho,d);
                % get predicted coordinates
                allpred = subtable.Cz{ho,1}(:,d);
                % index just coordinates in that holdout set
                pred = allpred(metadata(s).cvind(:,1)==ho);
                % calculate correlation
                output(ho,s,d) = corr(pred,utst);
                % calculate within-domain correlation
                animateoutput(ho,s,d) = corr(pred(1:5),utst(1:5));
                inanimateoutput(ho,s,d) = corr(pred(6:10),utst(6:10));
                if var(pred) == 0
                    output(ho,s,d) = 0;
                    animateoutput(ho,s,d) = 0;
                    inanimateoutput(ho,s,d) = 0;
                end     
            end
        end
    end
    result(r).output = output;
    result(r).animateoutput = animateoutput;
    result(r).inanimateoutput = inanimateoutput;

    % load perm results
    load([root,'/condor/derivatives/correlation/grOWL/performance/',ROI{r},'/perm_performance.mat']);
    
    % initialise output - !! participant x random seed x dimension (x fold, but we will
    % average over folds)
    output = zeros(size(subcode,2), 100, 3, 10);
    % initialise within-domain output
    animateoutput = zeros(size(subcode,2), 100, 3, 10);
    inanimateoutput = zeros(size(subcode,2), 100, 3, 10);

    % for each participant
    for s = 1:size(subcode,2)
        % get just the rows applying to that participant
        tmp = str2num(erase(subcode{s},'sub-'));
        subtable = Tallcv(Tallcv.subject==tmp,:);
        
        % for every random seed
        for seed = 1:100
            subsubtable = subtable(subtable.RandomSeed == seed,:);
            % for each holdout fold
            for ho = 1:10 
                % for each dimension
                for d = 1:3
                    % get target for that holdout set
                    utst = U(metadata(s).cvind(:,1)==ho,d);
                    % get predicted coordinates
                    allpred = subsubtable.Cz{ho,1}(:,d);
                    % index just coordinates in that holdout set
                    pred = allpred(metadata(s).cvind(:,1)==ho);
                    % calculate correlation
                    output(s,seed,d,ho) = corr(pred,utst);
                    % calculate within-domain correlation
                    animateoutput(s,seed,d,ho) = corr(pred(1:5),utst(1:5));
                    inanimateoutput(s,seed,d,ho) = corr(pred(6:10),utst(6:10));
                    if var(pred) == 0
                        output(s,seed,d,ho)= 0;
                        animateoutput(s,seed,d,ho) = 0;
                        inanimateoutput(s,seed,d,ho) = 0;
                    end
                end
            end
        end
    end
    permresult(r).output = output;
    permresult(r).animateoutput = animateoutput;
    permresult(r).inanimateoutput = inanimateoutput;
end

% bring output into workspace
whole_brain_grOWL = result(1).output;
left_ATL_grOWL = result(2).output;
left_ATL_ant_grOWL = result(3).output;
left_ATL_pos_grOWL = result(4).output;
right_ATL_grOWL = result(5).output;
right_ATL_ant_grOWL = result(6).output;
right_ATL_pos_grOWL = result(7).output;
perm_whole_brain_grOWL = mean(permresult(1).output,4);
perm_left_ATL_grOWL = mean(permresult(2).output,4);
perm_left_ATL_ant_grOWL = mean(permresult(3).output,4);
perm_left_ATL_pos_grOWL = mean(permresult(4).output,4);
perm_right_ATL_grOWL = mean(permresult(5).output,4);
perm_right_ATL_ant_grOWL = mean(permresult(6).output,4);
perm_right_ATL_pos_grOWL = mean(permresult(7).output,4);

whole_brain_grOWL_animate = result(1).animateoutput;
left_ATL_grOWL_animate = result(2).animateoutput;
left_ATL_ant_grOWL_animate = result(3).animateoutput;
left_ATL_pos_grOWL_animate = result(4).animateoutput;
right_ATL_grOWL_animate = result(5).animateoutput;
right_ATL_ant_grOWL_animate = result(6).animateoutput;
right_ATL_pos_grOWL_animate = result(7).animateoutput;
perm_whole_brain_grOWL_animate = mean(permresult(1).animateoutput,4);
perm_left_ATL_grOWL_animate = mean(permresult(2).animateoutput,4);
perm_left_ATL_ant_grOWL_animate = mean(permresult(3).animateoutput,4);
perm_left_ATL_pos_grOWL_animate = mean(permresult(4).animateoutput,4);
perm_right_ATL_grOWL_animate = mean(permresult(5).animateoutput,4);
perm_right_ATL_ant_grOWL_animate = mean(permresult(6).animateoutput,4);
perm_right_ATL_pos_grOWL_animate = mean(permresult(7).animateoutput,4);

whole_brain_grOWL_inanimate = result(1).inanimateoutput;
left_ATL_grOWL_inanimate = result(2).inanimateoutput;
left_ATL_ant_grOWL_inanimate = result(3).inanimateoutput;
left_ATL_pos_grOWL_inanimate = result(4).inanimateoutput;
right_ATL_grOWL_inanimate = result(5).inanimateoutput;
right_ATL_ant_grOWL_inanimate = result(6).inanimateoutput;
right_ATL_pos_grOWL_inanimate = result(7).inanimateoutput;
perm_whole_brain_grOWL_inanimate = mean(permresult(1).inanimateoutput,4);
perm_left_ATL_grOWL_inanimate = mean(permresult(2).inanimateoutput,4);
perm_left_ATL_ant_grOWL_inanimate = mean(permresult(3).inanimateoutput,4);
perm_left_ATL_pos_grOWL_inanimate = mean(permresult(4).inanimateoutput,4);
perm_right_ATL_grOWL_inanimate = mean(permresult(5).inanimateoutput,4);
perm_right_ATL_ant_grOWL_inanimate = mean(permresult(6).inanimateoutput,4);
perm_right_ATL_pos_grOWL_inanimate = mean(permresult(7).inanimateoutput,4);


% 2C. Statistics with Stelzer bootstrapping for the ATL (main analysis)

% initialise collated results
group_lin_LASSO_left = zeros(size(subcode,2),9);
group_lin_LASSO_right = zeros(size(subcode,2),9);
group_grOWL_left = zeros(size(subcode,2),9);
group_grOWL_right = zeros(size(subcode,2),9);

% collate
group_lin_LASSO_left(:,1) = mean(left_ATL_lin_LASSO_D1);
group_lin_LASSO_left(:,2) = mean(left_ATL_lin_LASSO_D2);
group_lin_LASSO_left(:,3) = mean(left_ATL_lin_LASSO_D3);
group_lin_LASSO_left(:,4) = mean(left_ATL_ant_lin_LASSO_D1);
group_lin_LASSO_left(:,5) = mean(left_ATL_ant_lin_LASSO_D2);
group_lin_LASSO_left(:,6) = mean(left_ATL_ant_lin_LASSO_D3);
group_lin_LASSO_left(:,7) = mean(left_ATL_pos_lin_LASSO_D1);
group_lin_LASSO_left(:,8) = mean(left_ATL_pos_lin_LASSO_D2);
group_lin_LASSO_left(:,9) = mean(left_ATL_pos_lin_LASSO_D3);

group_lin_LASSO_right(:,1) = mean(right_ATL_lin_LASSO_D1);
group_lin_LASSO_right(:,2) = mean(right_ATL_lin_LASSO_D2);
group_lin_LASSO_right(:,3) = mean(right_ATL_lin_LASSO_D3);
group_lin_LASSO_right(:,4) = mean(right_ATL_ant_lin_LASSO_D1);
group_lin_LASSO_right(:,5) = mean(right_ATL_ant_lin_LASSO_D2);
group_lin_LASSO_right(:,6) = mean(right_ATL_ant_lin_LASSO_D3);
group_lin_LASSO_right(:,7) = mean(right_ATL_pos_lin_LASSO_D1);
group_lin_LASSO_right(:,8) = mean(right_ATL_pos_lin_LASSO_D2);
group_lin_LASSO_right(:,9) = mean(right_ATL_pos_lin_LASSO_D3);

tmp = left_ATL_grOWL(:,:,1);
group_grOWL_left(:,1) = mean(tmp);
tmp = left_ATL_grOWL(:,:,2);
group_grOWL_left(:,2) = mean(tmp);
tmp = left_ATL_grOWL(:,:,3);
group_grOWL_left(:,3) = mean(tmp);
tmp = left_ATL_ant_grOWL(:,:,1);
group_grOWL_left(:,4) = mean(tmp);
tmp = left_ATL_ant_grOWL(:,:,2);
group_grOWL_left(:,5) = mean(tmp);
tmp = left_ATL_ant_grOWL(:,:,3);
group_grOWL_left(:,6) = mean(tmp);
tmp = left_ATL_pos_grOWL(:,:,1);
group_grOWL_left(:,7) = mean(tmp);
tmp = left_ATL_pos_grOWL(:,:,2);
group_grOWL_left(:,8) = mean(tmp);
tmp = left_ATL_pos_grOWL(:,:,3);
group_grOWL_left(:,9) = mean(tmp);

tmp = right_ATL_grOWL(:,:,1);
group_grOWL_right(:,1) = mean(tmp);
tmp = right_ATL_grOWL(:,:,2);
group_grOWL_right(:,2) = mean(tmp);
tmp = right_ATL_grOWL(:,:,3);
group_grOWL_right(:,3) = mean(tmp);
tmp = right_ATL_ant_grOWL(:,:,1);
group_grOWL_right(:,4) = mean(tmp);
tmp = right_ATL_ant_grOWL(:,:,2);
group_grOWL_right(:,5) = mean(tmp);
tmp = right_ATL_ant_grOWL(:,:,3);
group_grOWL_right(:,6) = mean(tmp);
tmp = right_ATL_pos_grOWL(:,:,1);
group_grOWL_right(:,7) = mean(tmp);
tmp = right_ATL_pos_grOWL(:,:,2);
group_grOWL_right(:,8) = mean(tmp);
tmp = right_ATL_pos_grOWL(:,:,3);
group_grOWL_right(:,9) = mean(tmp);

% construct permutation distributions by randomly sampling group averages
% (see Cox et al., 2024, Imaging Neuroscience for methods). Sorting makes
% it easier to caluclate percentile p values
bootstrapped_perm_lin_LASSO_left(:,1) = sort(mean(datasample(perm_left_ATL_lin_LASSO_D1,10000,2))');
bootstrapped_perm_lin_LASSO_left(:,2) = sort(mean(datasample(perm_left_ATL_lin_LASSO_D2,10000,2))');
bootstrapped_perm_lin_LASSO_left(:,3) = sort(mean(datasample(perm_left_ATL_lin_LASSO_D3,10000,2))');
bootstrapped_perm_lin_LASSO_left(:,4) = sort(mean(datasample(perm_left_ATL_ant_lin_LASSO_D1,10000,2))');
bootstrapped_perm_lin_LASSO_left(:,5) = sort(mean(datasample(perm_left_ATL_ant_lin_LASSO_D2,10000,2))');
bootstrapped_perm_lin_LASSO_left(:,6) = sort(mean(datasample(perm_left_ATL_ant_lin_LASSO_D3,10000,2))');
bootstrapped_perm_lin_LASSO_left(:,7) = sort(mean(datasample(perm_left_ATL_pos_lin_LASSO_D1,10000,2))');
bootstrapped_perm_lin_LASSO_left(:,8) = sort(mean(datasample(perm_left_ATL_pos_lin_LASSO_D2,10000,2))');
bootstrapped_perm_lin_LASSO_left(:,9) = sort(mean(datasample(perm_left_ATL_pos_lin_LASSO_D3,10000,2))');

bootstrapped_perm_lin_LASSO_right(:,1) = sort(mean(datasample(perm_right_ATL_lin_LASSO_D1,10000,2))');
bootstrapped_perm_lin_LASSO_right(:,2) = sort(mean(datasample(perm_right_ATL_lin_LASSO_D2,10000,2))');
bootstrapped_perm_lin_LASSO_right(:,3) = sort(mean(datasample(perm_right_ATL_lin_LASSO_D3,10000,2))');
bootstrapped_perm_lin_LASSO_right(:,4) = sort(mean(datasample(perm_right_ATL_ant_lin_LASSO_D1,10000,2))');
bootstrapped_perm_lin_LASSO_right(:,5) = sort(mean(datasample(perm_right_ATL_ant_lin_LASSO_D2,10000,2))');
bootstrapped_perm_lin_LASSO_right(:,6) = sort(mean(datasample(perm_right_ATL_ant_lin_LASSO_D3,10000,2))');
bootstrapped_perm_lin_LASSO_right(:,7) = sort(mean(datasample(perm_right_ATL_pos_lin_LASSO_D1,10000,2))');
bootstrapped_perm_lin_LASSO_right(:,8) = sort(mean(datasample(perm_right_ATL_pos_lin_LASSO_D2,10000,2))');
bootstrapped_perm_lin_LASSO_right(:,9) = sort(mean(datasample(perm_right_ATL_pos_lin_LASSO_D3,10000,2))');

tmp = perm_left_ATL_grOWL(:,:,1);
bootstrapped_perm_grOWL_left(:,1) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_grOWL(:,:,2);
bootstrapped_perm_grOWL_left(:,2) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_grOWL(:,:,3);
bootstrapped_perm_grOWL_left(:,3) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_ant_grOWL(:,:,1);
bootstrapped_perm_grOWL_left(:,4) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_ant_grOWL(:,:,2);
bootstrapped_perm_grOWL_left(:,5) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_ant_grOWL(:,:,3);
bootstrapped_perm_grOWL_left(:,6) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_pos_grOWL(:,:,1);
bootstrapped_perm_grOWL_left(:,7) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_pos_grOWL(:,:,2);
bootstrapped_perm_grOWL_left(:,8) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_pos_grOWL(:,:,3);
bootstrapped_perm_grOWL_left(:,9) = sort(mean(datasample(tmp,10000,2))');

tmp = perm_right_ATL_grOWL(:,:,1);
bootstrapped_perm_grOWL_right(:,1) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_grOWL(:,:,2);
bootstrapped_perm_grOWL_right(:,2) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_grOWL(:,:,3);
bootstrapped_perm_grOWL_right(:,3) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_ant_grOWL(:,:,1);
bootstrapped_perm_grOWL_right(:,4) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_ant_grOWL(:,:,2);
bootstrapped_perm_grOWL_right(:,5) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_ant_grOWL(:,:,3);
bootstrapped_perm_grOWL_right(:,6) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_pos_grOWL(:,:,1);
bootstrapped_perm_grOWL_right(:,7) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_pos_grOWL(:,:,2);
bootstrapped_perm_grOWL_right(:,8) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_pos_grOWL(:,:,3);
bootstrapped_perm_grOWL_right(:,9) = sort(mean(datasample(tmp,10000,2))');

correlation_p_val = zeros(4,9);
m = mean(group_lin_LASSO_left(:,1));
b = sum(bootstrapped_perm_lin_LASSO_left(:,1) > m);
correlation_p_val(1,1) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left(:,2));
b = sum(bootstrapped_perm_lin_LASSO_left(:,2) > m);
correlation_p_val(1,2) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left(:,3));
b = sum(bootstrapped_perm_lin_LASSO_left(:,3) > m);
correlation_p_val(1,3) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left(:,4));
b = sum(bootstrapped_perm_lin_LASSO_left(:,4) > m);
correlation_p_val(1,4) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left(:,5));
b = sum(bootstrapped_perm_lin_LASSO_left(:,5) > m);
correlation_p_val(1,5) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left(:,6));
b = sum(bootstrapped_perm_lin_LASSO_left(:,6) > m);
correlation_p_val(1,6) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left(:,7));
b = sum(bootstrapped_perm_lin_LASSO_left(:,7) > m);
correlation_p_val(1,7) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left(:,8));
b = sum(bootstrapped_perm_lin_LASSO_left(:,8) > m);
correlation_p_val(1,8) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left(:,9));
b = sum(bootstrapped_perm_lin_LASSO_left(:,9) > m);
correlation_p_val(1,9) = (b + 1)/(10000 + 1);

m = mean(group_lin_LASSO_right(:,1));
b = sum(bootstrapped_perm_lin_LASSO_right(:,1) > m);
correlation_p_val(2,1) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right(:,2));
b = sum(bootstrapped_perm_lin_LASSO_right(:,2) > m);
correlation_p_val(2,2) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right(:,3));
b = sum(bootstrapped_perm_lin_LASSO_right(:,3) > m);
correlation_p_val(2,3) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right(:,4));
b = sum(bootstrapped_perm_lin_LASSO_right(:,4) > m);
correlation_p_val(2,4) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right(:,5));
b = sum(bootstrapped_perm_lin_LASSO_right(:,5) > m);
correlation_p_val(2,5) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right(:,6));
b = sum(bootstrapped_perm_lin_LASSO_right(:,6) > m);
correlation_p_val(2,6) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right(:,7));
b = sum(bootstrapped_perm_lin_LASSO_right(:,7) > m);
correlation_p_val(2,7) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right(:,8));
b = sum(bootstrapped_perm_lin_LASSO_right(:,8) > m);
correlation_p_val(2,8) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right(:,9));
b = sum(bootstrapped_perm_lin_LASSO_right(:,9) > m);
correlation_p_val(2,9) = (b + 1)/(10000 + 1);

m = mean(group_grOWL_left(:,1));
b = sum(bootstrapped_perm_grOWL_left(:,1) > m);
correlation_p_val(3,1) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left(:,2));
b = sum(bootstrapped_perm_grOWL_left(:,2) > m);
correlation_p_val(3,2) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left(:,3));
b = sum(bootstrapped_perm_grOWL_left(:,3) > m);
correlation_p_val(3,3) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left(:,4));
b = sum(bootstrapped_perm_grOWL_left(:,4) > m);
correlation_p_val(3,4) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left(:,5));
b = sum(bootstrapped_perm_grOWL_left(:,5) > m);
correlation_p_val(3,5) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left(:,6));
b = sum(bootstrapped_perm_grOWL_left(:,6) > m);
correlation_p_val(3,6) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left(:,7));
b = sum(bootstrapped_perm_grOWL_left(:,7) > m);
correlation_p_val(3,7) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left(:,8));
b = sum(bootstrapped_perm_grOWL_left(:,8) > m);
correlation_p_val(3,8) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left(:,9));
b = sum(bootstrapped_perm_grOWL_left(:,9) > m);
correlation_p_val(3,9) = (b + 1)/(10000 + 1);

m = mean(group_grOWL_right(:,1));
b = sum(bootstrapped_perm_grOWL_right(:,1) > m);
correlation_p_val(4,1) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right(:,2));
b = sum(bootstrapped_perm_grOWL_right(:,2) > m);
correlation_p_val(4,2) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right(:,3));
b = sum(bootstrapped_perm_grOWL_right(:,3) > m);
correlation_p_val(4,3) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right(:,4));
b = sum(bootstrapped_perm_grOWL_right(:,4) > m);
correlation_p_val(4,4) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right(:,5));
b = sum(bootstrapped_perm_grOWL_right(:,5) > m);
correlation_p_val(4,5) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right(:,6));
b = sum(bootstrapped_perm_grOWL_right(:,6) > m);
correlation_p_val(4,6) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right(:,7));
b = sum(bootstrapped_perm_grOWL_right(:,7) > m);
correlation_p_val(4,7) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right(:,8));
b = sum(bootstrapped_perm_grOWL_right(:,8) > m);
correlation_p_val(4,8) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right(:,9));
b = sum(bootstrapped_perm_grOWL_right(:,9) > m);
correlation_p_val(4,9) = (b + 1)/(10000 + 1);

% control the false discovery rate
for i = 1:size(correlation_p_val,1)
    p = correlation_p_val(i,:);
    correlation_p_val(i,:) = mafdr(p,'BHFDR',1);
end

% test for differences between methods
correlation_difference_p_val = zeros(4,9);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_lin_LASSO_D1) - perm_left_ATL_grOWL(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_left(:,1) - group_grOWL_left(:,1));
b = sum(tmp > m);
correlation_difference_p_val(1,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_lin_LASSO_D2) - perm_left_ATL_grOWL(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_left(:,2) - group_grOWL_left(:,2));
b = sum(tmp > m);
correlation_difference_p_val(1,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_lin_LASSO_D3) - perm_left_ATL_grOWL(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_left(:,3) - group_grOWL_left(:,3));
b = sum(tmp > m);
correlation_difference_p_val(1,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_lin_LASSO_D1) - perm_left_ATL_ant_grOWL(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_left(:,4) - group_grOWL_left(:,4));
b = sum(tmp > m);
correlation_difference_p_val(1,4) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_lin_LASSO_D2) - perm_left_ATL_ant_grOWL(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_left(:,5) - group_grOWL_left(:,5));
b = sum(tmp > m);
correlation_difference_p_val(1,5) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_lin_LASSO_D3) - perm_left_ATL_ant_grOWL(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_left(:,6) - group_grOWL_left(:,6));
b = sum(tmp > m);
correlation_difference_p_val(1,6) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_lin_LASSO_D1) - perm_left_ATL_pos_grOWL(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_left(:,7) - group_grOWL_left(:,7));
b = sum(tmp > m);
correlation_difference_p_val(1,7) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_lin_LASSO_D2) - perm_left_ATL_pos_grOWL(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_left(:,8) - group_grOWL_left(:,8));
b = sum(tmp > m);
correlation_difference_p_val(1,8) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_lin_LASSO_D3) - perm_left_ATL_pos_grOWL(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_left(:,9) - group_grOWL_left(:,9));
b = sum(tmp > m);
correlation_difference_p_val(1,9) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_lin_LASSO_D1) - perm_right_ATL_grOWL(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_right(:,1) - group_grOWL_right(:,1));
b = sum(tmp > m);
correlation_difference_p_val(2,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_lin_LASSO_D2) - perm_right_ATL_grOWL(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_right(:,2) - group_grOWL_right(:,2));
b = sum(tmp > m);
correlation_difference_p_val(2,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_lin_LASSO_D3) - perm_right_ATL_grOWL(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_right(:,3) - group_grOWL_right(:,3));
b = sum(tmp > m);
correlation_difference_p_val(2,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_lin_LASSO_D1) - perm_right_ATL_ant_grOWL(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_right(:,4) - group_grOWL_right(:,4));
b = sum(tmp > m);
correlation_difference_p_val(2,4) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_lin_LASSO_D2) - perm_right_ATL_ant_grOWL(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_right(:,5) - group_grOWL_right(:,5));
b = sum(tmp > m);
correlation_difference_p_val(2,5) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_lin_LASSO_D3) - perm_right_ATL_ant_grOWL(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_right(:,6) - group_grOWL_right(:,6));
b = sum(tmp > m);
correlation_difference_p_val(2,6) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_lin_LASSO_D1) - perm_right_ATL_pos_grOWL(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_right(:,7) - group_grOWL_right(:,7));
b = sum(tmp > m);
correlation_difference_p_val(2,7) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_lin_LASSO_D2) - perm_right_ATL_pos_grOWL(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_right(:,8) - group_grOWL_right(:,8));
b = sum(tmp > m);
correlation_difference_p_val(2,8) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_lin_LASSO_D3) - perm_right_ATL_pos_grOWL(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_right(:,9) - group_grOWL_right(:,9));
b = sum(tmp > m);
correlation_difference_p_val(2,9) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_grOWL(:,:,1)) - perm_left_ATL_grOWL(:,:,1)),10000,2))');
m = mean(group_grOWL_left(:,1) - group_lin_LASSO_left(:,1));
b = sum(tmp > m);
correlation_difference_p_val(3,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_grOWL(:,:,2)) - perm_left_ATL_lin_LASSO_D2),10000,2))');
m = mean(group_grOWL_left(:,2) - group_lin_LASSO_left(:,2));
b = sum(tmp > m);
correlation_difference_p_val(3,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_grOWL(:,:,3)) - perm_left_ATL_lin_LASSO_D3),10000,2))');
m = mean(group_grOWL_left(:,3) - group_lin_LASSO_left(:,3));
b = sum(tmp > m);
correlation_difference_p_val(3,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_grOWL(:,:,1)) - perm_left_ATL_ant_lin_LASSO_D1),10000,2))');
m = mean(group_grOWL_left(:,4) - group_lin_LASSO_left(:,4));
b = sum(tmp > m);
correlation_difference_p_val(3,4) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_grOWL(:,:,2)) - perm_left_ATL_ant_lin_LASSO_D2),10000,2))');
m = mean(group_grOWL_left(:,5) - group_lin_LASSO_left(:,5));
b = sum(tmp > m);
correlation_difference_p_val(3,5) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_grOWL(:,:,3)) - perm_left_ATL_ant_lin_LASSO_D3),10000,2))');
m = mean(group_grOWL_left(:,6) - group_lin_LASSO_left(:,6));
b = sum(tmp > m);
correlation_difference_p_val(3,6) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_grOWL(:,:,1)) - perm_left_ATL_pos_lin_LASSO_D1),10000,2))');
m = mean(group_grOWL_left(:,7) - group_lin_LASSO_left(:,7));
b = sum(tmp > m);
correlation_difference_p_val(3,7) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_grOWL(:,:,2)) - perm_left_ATL_pos_lin_LASSO_D2),10000,2))');
m = mean(group_grOWL_left(:,8) - group_lin_LASSO_left(:,8));
b = sum(tmp > m);
correlation_difference_p_val(3,8) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_grOWL(:,:,3)) - perm_left_ATL_pos_lin_LASSO_D3),10000,2))');
m = mean(group_grOWL_left(:,9) - group_lin_LASSO_left(:,9));
b = sum(tmp > m);
correlation_difference_p_val(3,9) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_grOWL(:,:,1)) - perm_right_ATL_lin_LASSO_D1),10000,2))');
m = mean(group_grOWL_right(:,1) - group_lin_LASSO_right(:,1));
b = sum(tmp > m);
correlation_difference_p_val(4,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_grOWL(:,:,2)) - perm_right_ATL_lin_LASSO_D2),10000,2))');
m = mean(group_grOWL_right(:,2) - group_lin_LASSO_right(:,2));
b = sum(tmp > m);
correlation_difference_p_val(4,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_grOWL(:,:,3)) - perm_right_ATL_lin_LASSO_D3),10000,2))');
m = mean(group_grOWL_right(:,3) - group_lin_LASSO_right(:,3));
b = sum(tmp > m);
correlation_difference_p_val(4,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_grOWL(:,:,1)) - perm_right_ATL_ant_lin_LASSO_D1),10000,2))');
m = mean(group_grOWL_right(:,4) - group_lin_LASSO_right(:,4));
b = sum(tmp > m);
correlation_difference_p_val(4,4) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_grOWL(:,:,2)) - perm_right_ATL_ant_lin_LASSO_D2),10000,2))');
m = mean(group_grOWL_right(:,5) - group_lin_LASSO_right(:,5));
b = sum(tmp > m);
correlation_difference_p_val(4,5) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_grOWL(:,:,3)) - perm_right_ATL_ant_lin_LASSO_D3),10000,2))');
m = mean(group_grOWL_right(:,6) - group_lin_LASSO_right(:,6));
b = sum(tmp > m);
correlation_difference_p_val(4,6) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_grOWL(:,:,1)) - perm_right_ATL_pos_lin_LASSO_D1),10000,2))');
m = mean(group_grOWL_right(:,7) - group_lin_LASSO_right(:,7));
b = sum(tmp > m);
correlation_difference_p_val(4,7) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_grOWL(:,:,2)) - perm_right_ATL_pos_lin_LASSO_D2),10000,2))');
m = mean(group_grOWL_right(:,8) - group_lin_LASSO_right(:,8));
b = sum(tmp > m);
correlation_difference_p_val(4,8) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_grOWL(:,:,3)) - perm_right_ATL_pos_lin_LASSO_D3),10000,2))');
m = mean(group_grOWL_right(:,9) - group_lin_LASSO_right(:,9));
b = sum(tmp > m);
correlation_difference_p_val(4,9) = (b + 1)/(10000 + 1);

% control the false discovery rate
for i = 1:size(correlation_difference_p_val,1)
    p = correlation_difference_p_val(i,:);
    correlation_difference_p_val(i,:) = mafdr(p,'BHFDR',1);
end

% 2D. Plot results for the ATL (main analysis)

% plot
fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_lin_LASSO_left));
b.FaceColor = 'flat';
b.CData = colours_LASSO_left;
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_lin_LASSO_left)/sqrt(size(group_lin_LASSO_left,1)));
errorbar(1:size(group_lin_LASSO_left,2),mean(group_lin_LASSO_left,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(correlation_p_val,2)
    if (correlation_p_val(1,i) < 0.05 && correlation_difference_p_val(1,i) < 0.05)
       stars{i} = '*†';
    elseif correlation_p_val(1,i) < 0.05
        stars{i} = '*';
    elseif correlation_p_val(1,i) >= 0.05
       stars{i} = '';
    end
end
text((1:9)-0.1,mean(group_lin_LASSO_left)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/lin_LASSO_left.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_lin_LASSO_right));
b.FaceColor = 'flat';
b.CData = colours_LASSO_right;
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_lin_LASSO_right)/sqrt(size(group_lin_LASSO_right,1)));
errorbar(1:size(group_lin_LASSO_right,2),mean(group_lin_LASSO_right,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(correlation_p_val,2)
    if (correlation_p_val(2,i) < 0.05 && correlation_difference_p_val(2,i) < 0.05)
       stars{i} = '*†';
    elseif correlation_p_val(2,i) < 0.05
        stars{i} = '*';
    elseif correlation_p_val(2,i) >= 0.05
       stars{i} = '';
    end
end
text((1:9)-0.1,mean(group_lin_LASSO_right)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/lin_LASSO_right.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_grOWL_left));
b.FaceColor = 'flat';
b.CData = colours_grOWL_left;
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_grOWL_left)/sqrt(size(group_grOWL_left,1)));
errorbar(1:size(group_grOWL_left,2),mean(group_grOWL_left,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(correlation_p_val,2)
    if (correlation_p_val(3,i) < 0.05 && correlation_difference_p_val(3,i) < 0.05)
       stars{i} = '*†';
    elseif correlation_p_val(3,i) < 0.05
        stars{i} = '*';
    elseif correlation_p_val(3,i) >= 0.05
       stars{i} = '';
    end
end
text((1:9)-0.1,mean(group_grOWL_left)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/grOWL_left.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_grOWL_right));
b.FaceColor = 'flat';
b.CData = colours_grOWL_right;
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_grOWL_right)/sqrt(size(group_grOWL_right,1)));
errorbar(1:size(group_grOWL_right,2),mean(group_grOWL_right,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(correlation_p_val,2)
    if (correlation_p_val(4,i) < 0.05 && correlation_difference_p_val(4,i) < 0.05)
       stars{i} = '*†';
    elseif correlation_p_val(4,i) < 0.05
        stars{i} = '*';
    elseif correlation_p_val(4,i) >= 0.05
       stars{i} = '';
    end
end
text((1:9)-0.1,mean(group_grOWL_right)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/grOWL_right.png')
close(fig)

% 2E. Statistics with Stelzer bootstrapping for the ATL (within domain)

% initialise collated results
group_lin_LASSO_left_animate = zeros(size(subcode,2),9);
group_lin_LASSO_right_animate = zeros(size(subcode,2),9);
group_grOWL_left_animate = zeros(size(subcode,2),9);
group_grOWL_right_animate = zeros(size(subcode,2),9);
group_lin_LASSO_left_inanimate = zeros(size(subcode,2),9);
group_lin_LASSO_right_inanimate = zeros(size(subcode,2),9);
group_grOWL_left_inanimate = zeros(size(subcode,2),9);
group_grOWL_right_inanimate = zeros(size(subcode,2),9);

% collate
group_lin_LASSO_left_animate(:,1) = mean(left_ATL_lin_LASSO_D1_animate);
group_lin_LASSO_left_animate(:,2) = mean(left_ATL_lin_LASSO_D2_animate);
group_lin_LASSO_left_animate(:,3) = mean(left_ATL_lin_LASSO_D3_animate);
group_lin_LASSO_left_animate(:,4) = mean(left_ATL_ant_lin_LASSO_D1_animate);
group_lin_LASSO_left_animate(:,5) = mean(left_ATL_ant_lin_LASSO_D2_animate);
group_lin_LASSO_left_animate(:,6) = mean(left_ATL_ant_lin_LASSO_D3_animate);
group_lin_LASSO_left_animate(:,7) = mean(left_ATL_pos_lin_LASSO_D1_animate);
group_lin_LASSO_left_animate(:,8) = mean(left_ATL_pos_lin_LASSO_D2_animate);
group_lin_LASSO_left_animate(:,9) = mean(left_ATL_pos_lin_LASSO_D3_animate);

group_lin_LASSO_right_animate(:,1) = mean(right_ATL_lin_LASSO_D1_animate);
group_lin_LASSO_right_animate(:,2) = mean(right_ATL_lin_LASSO_D2_animate);
group_lin_LASSO_right_animate(:,3) = mean(right_ATL_lin_LASSO_D3_animate);
group_lin_LASSO_right_animate(:,4) = mean(right_ATL_ant_lin_LASSO_D1_animate);
group_lin_LASSO_right_animate(:,5) = mean(right_ATL_ant_lin_LASSO_D2_animate);
group_lin_LASSO_right_animate(:,6) = mean(right_ATL_ant_lin_LASSO_D3_animate);
group_lin_LASSO_right_animate(:,7) = mean(right_ATL_pos_lin_LASSO_D1_animate);
group_lin_LASSO_right_animate(:,8) = mean(right_ATL_pos_lin_LASSO_D2_animate);
group_lin_LASSO_right_animate(:,9) = mean(right_ATL_pos_lin_LASSO_D3_animate);

group_lin_LASSO_left_inanimate(:,1) = mean(left_ATL_lin_LASSO_D1_inanimate);
group_lin_LASSO_left_inanimate(:,2) = mean(left_ATL_lin_LASSO_D2_inanimate);
group_lin_LASSO_left_inanimate(:,3) = mean(left_ATL_lin_LASSO_D3_inanimate);
group_lin_LASSO_left_inanimate(:,4) = mean(left_ATL_ant_lin_LASSO_D1_inanimate);
group_lin_LASSO_left_inanimate(:,5) = mean(left_ATL_ant_lin_LASSO_D2_inanimate);
group_lin_LASSO_left_inanimate(:,6) = mean(left_ATL_ant_lin_LASSO_D3_inanimate);
group_lin_LASSO_left_inanimate(:,7) = mean(left_ATL_pos_lin_LASSO_D1_inanimate);
group_lin_LASSO_left_inanimate(:,8) = mean(left_ATL_pos_lin_LASSO_D2_inanimate);
group_lin_LASSO_left_inanimate(:,9) = mean(left_ATL_pos_lin_LASSO_D3_inanimate);

group_lin_LASSO_right_inanimate(:,1) = mean(right_ATL_lin_LASSO_D1_inanimate);
group_lin_LASSO_right_inanimate(:,2) = mean(right_ATL_lin_LASSO_D2_inanimate);
group_lin_LASSO_right_inanimate(:,3) = mean(right_ATL_lin_LASSO_D3_inanimate);
group_lin_LASSO_right_inanimate(:,4) = mean(right_ATL_ant_lin_LASSO_D1_inanimate);
group_lin_LASSO_right_inanimate(:,5) = mean(right_ATL_ant_lin_LASSO_D2_inanimate);
group_lin_LASSO_right_inanimate(:,6) = mean(right_ATL_ant_lin_LASSO_D3_inanimate);
group_lin_LASSO_right_inanimate(:,7) = mean(right_ATL_pos_lin_LASSO_D1_inanimate);
group_lin_LASSO_right_inanimate(:,8) = mean(right_ATL_pos_lin_LASSO_D2_inanimate);
group_lin_LASSO_right_inanimate(:,9) = mean(right_ATL_pos_lin_LASSO_D3_inanimate);

tmp = left_ATL_grOWL_animate(:,:,1);
group_grOWL_left_animate(:,1) = mean(tmp);
tmp = left_ATL_grOWL_animate(:,:,2);
group_grOWL_left_animate(:,2) = mean(tmp);
tmp = left_ATL_grOWL_animate(:,:,3);
group_grOWL_left_animate(:,3) = mean(tmp);
tmp = left_ATL_ant_grOWL_animate(:,:,1);
group_grOWL_left_animate(:,4) = mean(tmp);
tmp = left_ATL_ant_grOWL_animate(:,:,2);
group_grOWL_left_animate(:,5) = mean(tmp);
tmp = left_ATL_ant_grOWL_animate(:,:,3);
group_grOWL_left_animate(:,6) = mean(tmp);
tmp = left_ATL_pos_grOWL_animate(:,:,1);
group_grOWL_left_animate(:,7) = mean(tmp);
tmp = left_ATL_pos_grOWL_animate(:,:,2);
group_grOWL_left_animate(:,8) = mean(tmp);
tmp = left_ATL_pos_grOWL_animate(:,:,3);
group_grOWL_left_animate(:,9) = mean(tmp);

tmp = right_ATL_grOWL_animate(:,:,1);
group_grOWL_right_animate(:,1) = mean(tmp);
tmp = right_ATL_grOWL_animate(:,:,2);
group_grOWL_right_animate(:,2) = mean(tmp);
tmp = right_ATL_grOWL_animate(:,:,3);
group_grOWL_right_animate(:,3) = mean(tmp);
tmp = right_ATL_ant_grOWL_animate(:,:,1);
group_grOWL_right_animate(:,4) = mean(tmp);
tmp = right_ATL_ant_grOWL_animate(:,:,2);
group_grOWL_right_animate(:,5) = mean(tmp);
tmp = right_ATL_ant_grOWL_animate(:,:,3);
group_grOWL_right_animate(:,6) = mean(tmp);
tmp = right_ATL_pos_grOWL_animate(:,:,1);
group_grOWL_right_animate(:,7) = mean(tmp);
tmp = right_ATL_pos_grOWL_animate(:,:,2);
group_grOWL_right_animate(:,8) = mean(tmp);
tmp = right_ATL_pos_grOWL_animate(:,:,3);
group_grOWL_right_animate(:,9) = mean(tmp);

tmp = left_ATL_grOWL_inanimate(:,:,1);
group_grOWL_left_inanimate(:,1) = mean(tmp);
tmp = left_ATL_grOWL_inanimate(:,:,2);
group_grOWL_left_inanimate(:,2) = mean(tmp);
tmp = left_ATL_grOWL_inanimate(:,:,3);
group_grOWL_left_inanimate(:,3) = mean(tmp);
tmp = left_ATL_ant_grOWL_inanimate(:,:,1);
group_grOWL_left_inanimate(:,4) = mean(tmp);
tmp = left_ATL_ant_grOWL_inanimate(:,:,2);
group_grOWL_left_inanimate(:,5) = mean(tmp);
tmp = left_ATL_ant_grOWL_inanimate(:,:,3);
group_grOWL_left_inanimate(:,6) = mean(tmp);
tmp = left_ATL_pos_grOWL_inanimate(:,:,1);
group_grOWL_left_inanimate(:,7) = mean(tmp);
tmp = left_ATL_pos_grOWL_inanimate(:,:,2);
group_grOWL_left_inanimate(:,8) = mean(tmp);
tmp = left_ATL_pos_grOWL_inanimate(:,:,3);
group_grOWL_left_inanimate(:,9) = mean(tmp);

tmp = right_ATL_grOWL_inanimate(:,:,1);
group_grOWL_right_inanimate(:,1) = mean(tmp);
tmp = right_ATL_grOWL_inanimate(:,:,2);
group_grOWL_right_inanimate(:,2) = mean(tmp);
tmp = right_ATL_grOWL_inanimate(:,:,3);
group_grOWL_right_inanimate(:,3) = mean(tmp);
tmp = right_ATL_ant_grOWL_inanimate(:,:,1);
group_grOWL_right_inanimate(:,4) = mean(tmp);
tmp = right_ATL_ant_grOWL_inanimate(:,:,2);
group_grOWL_right_inanimate(:,5) = mean(tmp);
tmp = right_ATL_ant_grOWL_inanimate(:,:,3);
group_grOWL_right_inanimate(:,6) = mean(tmp);
tmp = right_ATL_pos_grOWL_inanimate(:,:,1);
group_grOWL_right_inanimate(:,7) = mean(tmp);
tmp = right_ATL_pos_grOWL_inanimate(:,:,2);
group_grOWL_right_inanimate(:,8) = mean(tmp);
tmp = right_ATL_pos_grOWL_inanimate(:,:,3);
group_grOWL_right_inanimate(:,9) = mean(tmp);

% construct permutation distributions by randomly sampling group averages
% (see Cox et al., 2024, Imaging Neuroscience for methods). Sorting makes
% it easier to caluclate percentile p values
bootstrapped_perm_lin_LASSO_left_animate(:,1) = sort(mean(datasample(perm_left_ATL_lin_LASSO_D1_animate,10000,2))');
bootstrapped_perm_lin_LASSO_left_animate(:,2) = sort(mean(datasample(perm_left_ATL_lin_LASSO_D2_animate,10000,2))');
bootstrapped_perm_lin_LASSO_left_animate(:,3) = sort(mean(datasample(perm_left_ATL_lin_LASSO_D3_animate,10000,2))');
bootstrapped_perm_lin_LASSO_left_animate(:,4) = sort(mean(datasample(perm_left_ATL_ant_lin_LASSO_D1_animate,10000,2))');
bootstrapped_perm_lin_LASSO_left_animate(:,5) = sort(mean(datasample(perm_left_ATL_ant_lin_LASSO_D2_animate,10000,2))');
bootstrapped_perm_lin_LASSO_left_animate(:,6) = sort(mean(datasample(perm_left_ATL_ant_lin_LASSO_D3_animate,10000,2))');
bootstrapped_perm_lin_LASSO_left_animate(:,7) = sort(mean(datasample(perm_left_ATL_pos_lin_LASSO_D1_animate,10000,2))');
bootstrapped_perm_lin_LASSO_left_animate(:,8) = sort(mean(datasample(perm_left_ATL_pos_lin_LASSO_D2_animate,10000,2))');
bootstrapped_perm_lin_LASSO_left_animate(:,9) = sort(mean(datasample(perm_left_ATL_pos_lin_LASSO_D3_animate,10000,2))');

bootstrapped_perm_lin_LASSO_right_animate(:,1) = sort(mean(datasample(perm_right_ATL_lin_LASSO_D1_animate,10000,2))');
bootstrapped_perm_lin_LASSO_right_animate(:,2) = sort(mean(datasample(perm_right_ATL_lin_LASSO_D2_animate,10000,2))');
bootstrapped_perm_lin_LASSO_right_animate(:,3) = sort(mean(datasample(perm_right_ATL_lin_LASSO_D3_animate,10000,2))');
bootstrapped_perm_lin_LASSO_right_animate(:,4) = sort(mean(datasample(perm_right_ATL_ant_lin_LASSO_D1_animate,10000,2))');
bootstrapped_perm_lin_LASSO_right_animate(:,5) = sort(mean(datasample(perm_right_ATL_ant_lin_LASSO_D2_animate,10000,2))');
bootstrapped_perm_lin_LASSO_right_animate(:,6) = sort(mean(datasample(perm_right_ATL_ant_lin_LASSO_D3_animate,10000,2))');
bootstrapped_perm_lin_LASSO_right_animate(:,7) = sort(mean(datasample(perm_right_ATL_pos_lin_LASSO_D1_animate,10000,2))');
bootstrapped_perm_lin_LASSO_right_animate(:,8) = sort(mean(datasample(perm_right_ATL_pos_lin_LASSO_D2_animate,10000,2))');
bootstrapped_perm_lin_LASSO_right_animate(:,9) = sort(mean(datasample(perm_right_ATL_pos_lin_LASSO_D3_animate,10000,2))');

bootstrapped_perm_lin_LASSO_left_inanimate(:,1) = sort(mean(datasample(perm_left_ATL_lin_LASSO_D1_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_left_inanimate(:,2) = sort(mean(datasample(perm_left_ATL_lin_LASSO_D2_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_left_inanimate(:,3) = sort(mean(datasample(perm_left_ATL_lin_LASSO_D3_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_left_inanimate(:,4) = sort(mean(datasample(perm_left_ATL_ant_lin_LASSO_D1_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_left_inanimate(:,5) = sort(mean(datasample(perm_left_ATL_ant_lin_LASSO_D2_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_left_inanimate(:,6) = sort(mean(datasample(perm_left_ATL_ant_lin_LASSO_D3_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_left_inanimate(:,7) = sort(mean(datasample(perm_left_ATL_pos_lin_LASSO_D1_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_left_inanimate(:,8) = sort(mean(datasample(perm_left_ATL_pos_lin_LASSO_D2_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_left_inanimate(:,9) = sort(mean(datasample(perm_left_ATL_pos_lin_LASSO_D3_inanimate,10000,2))');

bootstrapped_perm_lin_LASSO_right_inanimate(:,1) = sort(mean(datasample(perm_right_ATL_lin_LASSO_D1_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_right_inanimate(:,2) = sort(mean(datasample(perm_right_ATL_lin_LASSO_D2_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_right_inanimate(:,3) = sort(mean(datasample(perm_right_ATL_lin_LASSO_D3_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_right_inanimate(:,4) = sort(mean(datasample(perm_right_ATL_ant_lin_LASSO_D1_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_right_inanimate(:,5) = sort(mean(datasample(perm_right_ATL_ant_lin_LASSO_D2_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_right_inanimate(:,6) = sort(mean(datasample(perm_right_ATL_ant_lin_LASSO_D3_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_right_inanimate(:,7) = sort(mean(datasample(perm_right_ATL_pos_lin_LASSO_D1_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_right_inanimate(:,8) = sort(mean(datasample(perm_right_ATL_pos_lin_LASSO_D2_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_right_inanimate(:,9) = sort(mean(datasample(perm_right_ATL_pos_lin_LASSO_D3_inanimate,10000,2))');

tmp = perm_left_ATL_grOWL_animate(:,:,1);
bootstrapped_perm_grOWL_left_animate(:,1) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_grOWL_animate(:,:,2);
bootstrapped_perm_grOWL_left_animate(:,2) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_grOWL_animate(:,:,3);
bootstrapped_perm_grOWL_left_animate(:,3) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_ant_grOWL_animate(:,:,1);
bootstrapped_perm_grOWL_left_animate(:,4) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_ant_grOWL_animate(:,:,2);
bootstrapped_perm_grOWL_left_animate(:,5) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_ant_grOWL_animate(:,:,3);
bootstrapped_perm_grOWL_left_animate(:,6) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_pos_grOWL_animate(:,:,1);
bootstrapped_perm_grOWL_left_animate(:,7) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_pos_grOWL_animate(:,:,2);
bootstrapped_perm_grOWL_left_animate(:,8) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_pos_grOWL_animate(:,:,3);
bootstrapped_perm_grOWL_left_animate(:,9) = sort(mean(datasample(tmp,10000,2))');

tmp = perm_right_ATL_grOWL_animate(:,:,1);
bootstrapped_perm_grOWL_right_animate(:,1) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_grOWL_animate(:,:,2);
bootstrapped_perm_grOWL_right_animate(:,2) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_grOWL_animate(:,:,3);
bootstrapped_perm_grOWL_right_animate(:,3) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_ant_grOWL_animate(:,:,1);
bootstrapped_perm_grOWL_right_animate(:,4) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_ant_grOWL_animate(:,:,2);
bootstrapped_perm_grOWL_right_animate(:,5) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_ant_grOWL_animate(:,:,3);
bootstrapped_perm_grOWL_right_animate(:,6) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_pos_grOWL_animate(:,:,1);
bootstrapped_perm_grOWL_right_animate(:,7) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_pos_grOWL_animate(:,:,2);
bootstrapped_perm_grOWL_right_animate(:,8) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_pos_grOWL_animate(:,:,3);
bootstrapped_perm_grOWL_right_animate(:,9) = sort(mean(datasample(tmp,10000,2))');

tmp = perm_left_ATL_grOWL_inanimate(:,:,1);
bootstrapped_perm_grOWL_left_inanimate(:,1) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_grOWL_inanimate(:,:,2);
bootstrapped_perm_grOWL_left_inanimate(:,2) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_grOWL_inanimate(:,:,3);
bootstrapped_perm_grOWL_left_inanimate(:,3) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_ant_grOWL_inanimate(:,:,1);
bootstrapped_perm_grOWL_left_inanimate(:,4) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_ant_grOWL_inanimate(:,:,2);
bootstrapped_perm_grOWL_left_inanimate(:,5) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_ant_grOWL_inanimate(:,:,3);
bootstrapped_perm_grOWL_left_inanimate(:,6) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_pos_grOWL_inanimate(:,:,1);
bootstrapped_perm_grOWL_left_inanimate(:,7) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_pos_grOWL_inanimate(:,:,2);
bootstrapped_perm_grOWL_left_inanimate(:,8) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_left_ATL_pos_grOWL_inanimate(:,:,3);
bootstrapped_perm_grOWL_left_inanimate(:,9) = sort(mean(datasample(tmp,10000,2))');

tmp = perm_right_ATL_grOWL_inanimate(:,:,1);
bootstrapped_perm_grOWL_right_inanimate(:,1) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_grOWL_inanimate(:,:,2);
bootstrapped_perm_grOWL_right_inanimate(:,2) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_grOWL_inanimate(:,:,3);
bootstrapped_perm_grOWL_right_inanimate(:,3) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_ant_grOWL_inanimate(:,:,1);
bootstrapped_perm_grOWL_right_inanimate(:,4) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_ant_grOWL_inanimate(:,:,2);
bootstrapped_perm_grOWL_right_inanimate(:,5) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_ant_grOWL_inanimate(:,:,3);
bootstrapped_perm_grOWL_right_inanimate(:,6) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_pos_grOWL_inanimate(:,:,1);
bootstrapped_perm_grOWL_right_inanimate(:,7) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_pos_grOWL_inanimate(:,:,2);
bootstrapped_perm_grOWL_right_inanimate(:,8) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_right_ATL_pos_grOWL_inanimate(:,:,3);
bootstrapped_perm_grOWL_right_inanimate(:,9) = sort(mean(datasample(tmp,10000,2))');

animate_correlation_p_val = zeros(4,9);
m = mean(group_lin_LASSO_left_animate(:,1));
b = sum(bootstrapped_perm_lin_LASSO_left_animate(:,1) > m);
animate_correlation_p_val(1,1) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_animate(:,2));
b = sum(bootstrapped_perm_lin_LASSO_left_animate(:,2) > m);
animate_correlation_p_val(1,2) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_animate(:,3));
b = sum(bootstrapped_perm_lin_LASSO_left_animate(:,3) > m);
animate_correlation_p_val(1,3) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_animate(:,4));
b = sum(bootstrapped_perm_lin_LASSO_left_animate(:,4) > m);
animate_correlation_p_val(1,4) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_animate(:,5));
b = sum(bootstrapped_perm_lin_LASSO_left_animate(:,5) > m);
animate_correlation_p_val(1,5) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_animate(:,6));
b = sum(bootstrapped_perm_lin_LASSO_left_animate(:,6) > m);
animate_correlation_p_val(1,6) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_animate(:,7));
b = sum(bootstrapped_perm_lin_LASSO_left_animate(:,7) > m);
animate_correlation_p_val(1,7) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_animate(:,8));
b = sum(bootstrapped_perm_lin_LASSO_left_animate(:,8) > m);
animate_correlation_p_val(1,8) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_animate(:,9));
b = sum(bootstrapped_perm_lin_LASSO_left_animate(:,9) > m);
animate_correlation_p_val(1,9) = (b + 1)/(10000 + 1);

m = mean(group_lin_LASSO_right_animate(:,1));
b = sum(bootstrapped_perm_lin_LASSO_right_animate(:,1) > m);
animate_correlation_p_val(2,1) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_animate(:,2));
b = sum(bootstrapped_perm_lin_LASSO_right_animate(:,2) > m);
animate_correlation_p_val(2,2) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_animate(:,3));
b = sum(bootstrapped_perm_lin_LASSO_right_animate(:,3) > m);
animate_correlation_p_val(2,3) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_animate(:,4));
b = sum(bootstrapped_perm_lin_LASSO_right_animate(:,4) > m);
animate_correlation_p_val(2,4) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_animate(:,5));
b = sum(bootstrapped_perm_lin_LASSO_right_animate(:,5) > m);
animate_correlation_p_val(2,5) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_animate(:,6));
b = sum(bootstrapped_perm_lin_LASSO_right_animate(:,6) > m);
animate_correlation_p_val(2,6) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_animate(:,7));
b = sum(bootstrapped_perm_lin_LASSO_right_animate(:,7) > m);
animate_correlation_p_val(2,7) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_animate(:,8));
b = sum(bootstrapped_perm_lin_LASSO_right_animate(:,8) > m);
animate_correlation_p_val(2,8) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_animate(:,9));
b = sum(bootstrapped_perm_lin_LASSO_right_animate(:,9) > m);
animate_correlation_p_val(2,9) = (b + 1)/(10000 + 1);

m = mean(group_grOWL_left_animate(:,1));
b = sum(bootstrapped_perm_grOWL_left_animate(:,1) > m);
animate_correlation_p_val(3,1) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_animate(:,2));
b = sum(bootstrapped_perm_grOWL_left_animate(:,2) > m);
animate_correlation_p_val(3,2) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_animate(:,3));
b = sum(bootstrapped_perm_grOWL_left_animate(:,3) > m);
animate_correlation_p_val(3,3) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_animate(:,4));
b = sum(bootstrapped_perm_grOWL_left_animate(:,4) > m);
animate_correlation_p_val(3,4) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_animate(:,5));
b = sum(bootstrapped_perm_grOWL_left_animate(:,5) > m);
animate_correlation_p_val(3,5) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_animate(:,6));
b = sum(bootstrapped_perm_grOWL_left_animate(:,6) > m);
animate_correlation_p_val(3,6) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_animate(:,7));
b = sum(bootstrapped_perm_grOWL_left_animate(:,7) > m);
animate_correlation_p_val(3,7) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_animate(:,8));
b = sum(bootstrapped_perm_grOWL_left_animate(:,8) > m);
animate_correlation_p_val(3,8) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_animate(:,9));
b = sum(bootstrapped_perm_grOWL_left_animate(:,9) > m);
animate_correlation_p_val(3,9) = (b + 1)/(10000 + 1);

m = mean(group_grOWL_right_animate(:,1));
b = sum(bootstrapped_perm_grOWL_right_animate(:,1) > m);
animate_correlation_p_val(4,1) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_animate(:,2));
b = sum(bootstrapped_perm_grOWL_right_animate(:,2) > m);
animate_correlation_p_val(4,2) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_animate(:,3));
b = sum(bootstrapped_perm_grOWL_right_animate(:,3) > m);
animate_correlation_p_val(4,3) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_animate(:,4));
b = sum(bootstrapped_perm_grOWL_right_animate(:,4) > m);
animate_correlation_p_val(4,4) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_animate(:,5));
b = sum(bootstrapped_perm_grOWL_right_animate(:,5) > m);
animate_correlation_p_val(4,5) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_animate(:,6));
b = sum(bootstrapped_perm_grOWL_right_animate(:,6) > m);
animate_correlation_p_val(4,6) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_animate(:,7));
b = sum(bootstrapped_perm_grOWL_right_animate(:,7) > m);
animate_correlation_p_val(4,7) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_animate(:,8));
b = sum(bootstrapped_perm_grOWL_right_animate(:,8) > m);
animate_correlation_p_val(4,8) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_animate(:,9));
b = sum(bootstrapped_perm_grOWL_right_animate(:,9) > m);
animate_correlation_p_val(4,9) = (b + 1)/(10000 + 1);

inanimate_correlation_p_val = zeros(4,9);
m = mean(group_lin_LASSO_left_inanimate(:,1));
b = sum(bootstrapped_perm_lin_LASSO_left_inanimate(:,1) > m);
inanimate_correlation_p_val(1,1) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_inanimate(:,2));
b = sum(bootstrapped_perm_lin_LASSO_left_inanimate(:,2) > m);
inanimate_correlation_p_val(1,2) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_inanimate(:,3));
b = sum(bootstrapped_perm_lin_LASSO_left_inanimate(:,3) > m);
inanimate_correlation_p_val(1,3) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_inanimate(:,4));
b = sum(bootstrapped_perm_lin_LASSO_left_inanimate(:,4) > m);
inanimate_correlation_p_val(1,4) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_inanimate(:,5));
b = sum(bootstrapped_perm_lin_LASSO_left_inanimate(:,5) > m);
inanimate_correlation_p_val(1,5) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_inanimate(:,6));
b = sum(bootstrapped_perm_lin_LASSO_left_inanimate(:,6) > m);
inanimate_correlation_p_val(1,6) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_inanimate(:,7));
b = sum(bootstrapped_perm_lin_LASSO_left_inanimate(:,7) > m);
inanimate_correlation_p_val(1,7) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_inanimate(:,8));
b = sum(bootstrapped_perm_lin_LASSO_left_inanimate(:,8) > m);
inanimate_correlation_p_val(1,8) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_left_inanimate(:,9));
b = sum(bootstrapped_perm_lin_LASSO_left_inanimate(:,9) > m);
inanimate_correlation_p_val(1,9) = (b + 1)/(10000 + 1);

m = mean(group_lin_LASSO_right_inanimate(:,1));
b = sum(bootstrapped_perm_lin_LASSO_right_inanimate(:,1) > m);
inanimate_correlation_p_val(2,1) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_inanimate(:,2));
b = sum(bootstrapped_perm_lin_LASSO_right_inanimate(:,2) > m);
inanimate_correlation_p_val(2,2) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_inanimate(:,3));
b = sum(bootstrapped_perm_lin_LASSO_right_inanimate(:,3) > m);
inanimate_correlation_p_val(2,3) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_inanimate(:,4));
b = sum(bootstrapped_perm_lin_LASSO_right_inanimate(:,4) > m);
inanimate_correlation_p_val(2,4) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_inanimate(:,5));
b = sum(bootstrapped_perm_lin_LASSO_right_inanimate(:,5) > m);
inanimate_correlation_p_val(2,5) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_inanimate(:,6));
b = sum(bootstrapped_perm_lin_LASSO_right_inanimate(:,6) > m);
inanimate_correlation_p_val(2,6) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_inanimate(:,7));
b = sum(bootstrapped_perm_lin_LASSO_right_inanimate(:,7) > m);
inanimate_correlation_p_val(2,7) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_inanimate(:,8));
b = sum(bootstrapped_perm_lin_LASSO_right_inanimate(:,8) > m);
inanimate_correlation_p_val(2,8) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_right_inanimate(:,9));
b = sum(bootstrapped_perm_lin_LASSO_right_inanimate(:,9) > m);
inanimate_correlation_p_val(2,9) = (b + 1)/(10000 + 1);

m = mean(group_grOWL_left_inanimate(:,1));
b = sum(bootstrapped_perm_grOWL_left_inanimate(:,1) > m);
inanimate_correlation_p_val(3,1) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_inanimate(:,2));
b = sum(bootstrapped_perm_grOWL_left_inanimate(:,2) > m);
inanimate_correlation_p_val(3,2) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_inanimate(:,3));
b = sum(bootstrapped_perm_grOWL_left_inanimate(:,3) > m);
inanimate_correlation_p_val(3,3) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_inanimate(:,4));
b = sum(bootstrapped_perm_grOWL_left_inanimate(:,4) > m);
inanimate_correlation_p_val(3,4) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_inanimate(:,5));
b = sum(bootstrapped_perm_grOWL_left_inanimate(:,5) > m);
inanimate_correlation_p_val(3,5) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_inanimate(:,6));
b = sum(bootstrapped_perm_grOWL_left_inanimate(:,6) > m);
inanimate_correlation_p_val(3,6) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_inanimate(:,7));
b = sum(bootstrapped_perm_grOWL_left_inanimate(:,7) > m);
inanimate_correlation_p_val(3,7) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_inanimate(:,8));
b = sum(bootstrapped_perm_grOWL_left_inanimate(:,8) > m);
inanimate_correlation_p_val(3,8) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_left_inanimate(:,9));
b = sum(bootstrapped_perm_grOWL_left_inanimate(:,9) > m);
inanimate_correlation_p_val(3,9) = (b + 1)/(10000 + 1);

m = mean(group_grOWL_right_inanimate(:,1));
b = sum(bootstrapped_perm_grOWL_right_inanimate(:,1) > m);
inanimate_correlation_p_val(4,1) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_inanimate(:,2));
b = sum(bootstrapped_perm_grOWL_right_inanimate(:,2) > m);
inanimate_correlation_p_val(4,2) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_inanimate(:,3));
b = sum(bootstrapped_perm_grOWL_right_inanimate(:,3) > m);
inanimate_correlation_p_val(4,3) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_inanimate(:,4));
b = sum(bootstrapped_perm_grOWL_right_inanimate(:,4) > m);
inanimate_correlation_p_val(4,4) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_inanimate(:,5));
b = sum(bootstrapped_perm_grOWL_right_inanimate(:,5) > m);
inanimate_correlation_p_val(4,5) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_inanimate(:,6));
b = sum(bootstrapped_perm_grOWL_right_inanimate(:,6) > m);
inanimate_correlation_p_val(4,6) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_inanimate(:,7));
b = sum(bootstrapped_perm_grOWL_right_inanimate(:,7) > m);
inanimate_correlation_p_val(4,7) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_inanimate(:,8));
b = sum(bootstrapped_perm_grOWL_right_inanimate(:,8) > m);
inanimate_correlation_p_val(4,8) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_right_inanimate(:,9));
b = sum(bootstrapped_perm_grOWL_right_inanimate(:,9) > m);
inanimate_correlation_p_val(4,9) = (b + 1)/(10000 + 1);

% test for differences between methods
animate_correlation_difference_p_val = zeros(4,9);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_lin_LASSO_D1_animate) - perm_left_ATL_grOWL_animate(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_left_animate(:,1) - group_grOWL_left_animate(:,1));
b = sum(tmp > m);
animate_correlation_difference_p_val(1,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_lin_LASSO_D2_animate) - perm_left_ATL_grOWL_animate(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_left_animate(:,2) - group_grOWL_left_animate(:,2));
b = sum(tmp > m);
animate_correlation_difference_p_val(1,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_lin_LASSO_D3_animate) - perm_left_ATL_grOWL_animate(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_left_animate(:,3) - group_grOWL_left_animate(:,3));
b = sum(tmp > m);
animate_correlation_difference_p_val(1,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_lin_LASSO_D1_animate) - perm_left_ATL_ant_grOWL_animate(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_left_animate(:,4) - group_grOWL_left_animate(:,4));
b = sum(tmp > m);
animate_correlation_difference_p_val(1,4) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_lin_LASSO_D2_animate) - perm_left_ATL_ant_grOWL_animate(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_left_animate(:,5) - group_grOWL_left_animate(:,5));
b = sum(tmp > m);
animate_correlation_difference_p_val(1,5) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_lin_LASSO_D3_animate) - perm_left_ATL_ant_grOWL_animate(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_left_animate(:,6) - group_grOWL_left_animate(:,6));
b = sum(tmp > m);
animate_correlation_difference_p_val(1,6) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_lin_LASSO_D1_animate) - perm_left_ATL_pos_grOWL_animate(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_left_animate(:,7) - group_grOWL_left_animate(:,7));
b = sum(tmp > m);
animate_correlation_difference_p_val(1,7) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_lin_LASSO_D2_animate) - perm_left_ATL_pos_grOWL_animate(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_left_animate(:,8) - group_grOWL_left_animate(:,8));
b = sum(tmp > m);
animate_correlation_difference_p_val(1,8) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_lin_LASSO_D3_animate) - perm_left_ATL_pos_grOWL_animate(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_left_animate(:,9) - group_grOWL_left_animate(:,9));
b = sum(tmp > m);
animate_correlation_difference_p_val(1,9) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_lin_LASSO_D1_animate) - perm_right_ATL_grOWL_animate(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_right_animate(:,1) - group_grOWL_right_animate(:,1));
b = sum(tmp > m);
animate_correlation_difference_p_val(2,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_lin_LASSO_D2_animate) - perm_right_ATL_grOWL_animate(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_right_animate(:,2) - group_grOWL_right_animate(:,2));
b = sum(tmp > m);
animate_correlation_difference_p_val(2,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_lin_LASSO_D3_animate) - perm_right_ATL_grOWL_animate(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_right_animate(:,3) - group_grOWL_right_animate(:,3));
b = sum(tmp > m);
animate_correlation_difference_p_val(2,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_lin_LASSO_D1_animate) - perm_right_ATL_ant_grOWL_animate(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_right_animate(:,4) - group_grOWL_right_animate(:,4));
b = sum(tmp > m);
animate_correlation_difference_p_val(2,4) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_lin_LASSO_D2_animate) - perm_right_ATL_ant_grOWL_animate(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_right_animate(:,5) - group_grOWL_right_animate(:,5));
b = sum(tmp > m);
animate_correlation_difference_p_val(2,5) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_lin_LASSO_D3_animate) - perm_right_ATL_ant_grOWL_animate(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_right_animate(:,6) - group_grOWL_right_animate(:,6));
b = sum(tmp > m);
animate_correlation_difference_p_val(2,6) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_lin_LASSO_D1_animate) - perm_right_ATL_pos_grOWL_animate(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_right_animate(:,7) - group_grOWL_right_animate(:,7));
b = sum(tmp > m);
animate_correlation_difference_p_val(2,7) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_lin_LASSO_D2_animate) - perm_right_ATL_pos_grOWL_animate(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_right_animate(:,8) - group_grOWL_right_animate(:,8));
b = sum(tmp > m);
animate_correlation_difference_p_val(2,8) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_lin_LASSO_D3_animate) - perm_right_ATL_pos_grOWL_animate(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_right_animate(:,9) - group_grOWL_right_animate(:,9));
b = sum(tmp > m);
animate_correlation_difference_p_val(2,9) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_grOWL_animate(:,:,1)) - perm_left_ATL_grOWL_animate(:,:,1)),10000,2))');
m = mean(group_grOWL_left_animate(:,1) - group_lin_LASSO_left_animate(:,1));
b = sum(tmp > m);
animate_correlation_difference_p_val(3,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_grOWL_animate(:,:,2)) - perm_left_ATL_lin_LASSO_D2_animate),10000,2))');
m = mean(group_grOWL_left_animate(:,2) - group_lin_LASSO_left_animate(:,2));
b = sum(tmp > m);
animate_correlation_difference_p_val(3,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_grOWL_animate(:,:,3)) - perm_left_ATL_lin_LASSO_D3_animate),10000,2))');
m = mean(group_grOWL_left_animate(:,3) - group_lin_LASSO_left_animate(:,3));
b = sum(tmp > m);
animate_correlation_difference_p_val(3,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_grOWL_animate(:,:,1)) - perm_left_ATL_ant_lin_LASSO_D1_animate),10000,2))');
m = mean(group_grOWL_left_animate(:,4) - group_lin_LASSO_left_animate(:,4));
b = sum(tmp > m);
animate_correlation_difference_p_val(3,4) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_grOWL_animate(:,:,2)) - perm_left_ATL_ant_lin_LASSO_D2_animate),10000,2))');
m = mean(group_grOWL_left_animate(:,5) - group_lin_LASSO_left_animate(:,5));
b = sum(tmp > m);
animate_correlation_difference_p_val(3,5) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_grOWL_animate(:,:,3)) - perm_left_ATL_ant_lin_LASSO_D3_animate),10000,2))');
m = mean(group_grOWL_left_animate(:,6) - group_lin_LASSO_left_animate(:,6));
b = sum(tmp > m);
animate_correlation_difference_p_val(3,6) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_grOWL_animate(:,:,1)) - perm_left_ATL_pos_lin_LASSO_D1_animate),10000,2))');
m = mean(group_grOWL_left_animate(:,7) - group_lin_LASSO_left_animate(:,7));
b = sum(tmp > m);
animate_correlation_difference_p_val(3,7) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_grOWL_animate(:,:,2)) - perm_left_ATL_pos_lin_LASSO_D2_animate),10000,2))');
m = mean(group_grOWL_left_animate(:,8) - group_lin_LASSO_left_animate(:,8));
b = sum(tmp > m);
animate_correlation_difference_p_val(3,8) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_grOWL_animate(:,:,3)) - perm_left_ATL_pos_lin_LASSO_D3_animate),10000,2))');
m = mean(group_grOWL_left_animate(:,9) - group_lin_LASSO_left_animate(:,9));
b = sum(tmp > m);
animate_correlation_difference_p_val(3,9) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_grOWL_animate(:,:,1)) - perm_right_ATL_lin_LASSO_D1_animate),10000,2))');
m = mean(group_grOWL_right_animate(:,1) - group_lin_LASSO_right_animate(:,1));
b = sum(tmp > m);
animate_correlation_difference_p_val(4,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_grOWL_animate(:,:,2)) - perm_right_ATL_lin_LASSO_D2_animate),10000,2))');
m = mean(group_grOWL_right_animate(:,2) - group_lin_LASSO_right_animate(:,2));
b = sum(tmp > m);
animate_correlation_difference_p_val(4,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_grOWL_animate(:,:,3)) - perm_right_ATL_lin_LASSO_D3_animate),10000,2))');
m = mean(group_grOWL_right_animate(:,3) - group_lin_LASSO_right_animate(:,3));
b = sum(tmp > m);
animate_correlation_difference_p_val(4,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_grOWL_animate(:,:,1)) - perm_right_ATL_ant_lin_LASSO_D1_animate),10000,2))');
m = mean(group_grOWL_right_animate(:,4) - group_lin_LASSO_right_animate(:,4));
b = sum(tmp > m);
animate_correlation_difference_p_val(4,4) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_grOWL_animate(:,:,2)) - perm_right_ATL_ant_lin_LASSO_D2_animate),10000,2))');
m = mean(group_grOWL_right_animate(:,5) - group_lin_LASSO_right_animate(:,5));
b = sum(tmp > m);
animate_correlation_difference_p_val(4,5) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_grOWL_animate(:,:,3)) - perm_right_ATL_ant_lin_LASSO_D3_animate),10000,2))');
m = mean(group_grOWL_right_animate(:,6) - group_lin_LASSO_right_animate(:,6));
b = sum(tmp > m);
animate_correlation_difference_p_val(4,6) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_grOWL_animate(:,:,1)) - perm_right_ATL_pos_lin_LASSO_D1_animate),10000,2))');
m = mean(group_grOWL_right_animate(:,7) - group_lin_LASSO_right_animate(:,7));
b = sum(tmp > m);
animate_correlation_difference_p_val(4,7) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_grOWL_animate(:,:,2)) - perm_right_ATL_pos_lin_LASSO_D2_animate),10000,2))');
m = mean(group_grOWL_right_animate(:,8) - group_lin_LASSO_right_animate(:,8));
b = sum(tmp > m);
animate_correlation_difference_p_val(4,8) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_grOWL_animate(:,:,3)) - perm_right_ATL_pos_lin_LASSO_D3_animate),10000,2))');
m = mean(group_grOWL_right_animate(:,9) - group_lin_LASSO_right_animate(:,9));
b = sum(tmp > m);
animate_correlation_difference_p_val(4,9) = (b + 1)/(10000 + 1);


% control the false discovery rate
for i = 1:size(animate_correlation_difference_p_val,1)
    p = animate_correlation_difference_p_val(i,:);
    animate_correlation_difference_p_val(i,:) = mafdr(p,'BHFDR',1);
end

inanimate_correlation_difference_p_val = zeros(4,9);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_lin_LASSO_D1_inanimate) - perm_left_ATL_grOWL_inanimate(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_left_inanimate(:,1) - group_grOWL_left_inanimate(:,1));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(1,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_lin_LASSO_D2_inanimate) - perm_left_ATL_grOWL_inanimate(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_left_inanimate(:,2) - group_grOWL_left_inanimate(:,2));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(1,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_lin_LASSO_D3_inanimate) - perm_left_ATL_grOWL_inanimate(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_left_inanimate(:,3) - group_grOWL_left_inanimate(:,3));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(1,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_lin_LASSO_D1_inanimate) - perm_left_ATL_ant_grOWL_inanimate(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_left_inanimate(:,4) - group_grOWL_left_inanimate(:,4));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(1,4) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_lin_LASSO_D2_inanimate) - perm_left_ATL_ant_grOWL_inanimate(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_left_inanimate(:,5) - group_grOWL_left_inanimate(:,5));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(1,5) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_lin_LASSO_D3_inanimate) - perm_left_ATL_ant_grOWL_inanimate(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_left_inanimate(:,6) - group_grOWL_left_inanimate(:,6));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(1,6) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_lin_LASSO_D1_inanimate) - perm_left_ATL_pos_grOWL_inanimate(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_left_inanimate(:,7) - group_grOWL_left_inanimate(:,7));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(1,7) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_lin_LASSO_D2_inanimate) - perm_left_ATL_pos_grOWL_inanimate(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_left_inanimate(:,8) - group_grOWL_left_inanimate(:,8));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(1,8) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_lin_LASSO_D3_inanimate) - perm_left_ATL_pos_grOWL_inanimate(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_left_inanimate(:,9) - group_grOWL_left_inanimate(:,9));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(1,9) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_lin_LASSO_D1_inanimate) - perm_right_ATL_grOWL_inanimate(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_right_inanimate(:,1) - group_grOWL_right_inanimate(:,1));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(2,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_lin_LASSO_D2_inanimate) - perm_right_ATL_grOWL_inanimate(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_right_inanimate(:,2) - group_grOWL_right_inanimate(:,2));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(2,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_lin_LASSO_D3_inanimate) - perm_right_ATL_grOWL_inanimate(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_right_inanimate(:,3) - group_grOWL_right_inanimate(:,3));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(2,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_lin_LASSO_D1_inanimate) - perm_right_ATL_ant_grOWL_inanimate(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_right_inanimate(:,4) - group_grOWL_right_inanimate(:,4));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(2,4) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_lin_LASSO_D2_inanimate) - perm_right_ATL_ant_grOWL_inanimate(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_right_inanimate(:,5) - group_grOWL_right_inanimate(:,5));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(2,5) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_lin_LASSO_D3_inanimate) - perm_right_ATL_ant_grOWL_inanimate(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_right_inanimate(:,6) - group_grOWL_right_inanimate(:,6));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(2,6) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_lin_LASSO_D1_inanimate) - perm_right_ATL_pos_grOWL_inanimate(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_right_inanimate(:,7) - group_grOWL_right_inanimate(:,7));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(2,7) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_lin_LASSO_D2_inanimate) - perm_right_ATL_pos_grOWL_inanimate(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_right_inanimate(:,8) - group_grOWL_right_inanimate(:,8));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(2,8) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_lin_LASSO_D3_inanimate) - perm_right_ATL_pos_grOWL_inanimate(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_right_inanimate(:,9) - group_grOWL_right_inanimate(:,9));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(2,9) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_grOWL_inanimate(:,:,1)) - perm_left_ATL_grOWL_inanimate(:,:,1)),10000,2))');
m = mean(group_grOWL_left_inanimate(:,1) - group_lin_LASSO_left_inanimate(:,1));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(3,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_grOWL_inanimate(:,:,2)) - perm_left_ATL_lin_LASSO_D2_inanimate),10000,2))');
m = mean(group_grOWL_left_inanimate(:,2) - group_lin_LASSO_left_inanimate(:,2));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(3,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_grOWL_inanimate(:,:,3)) - perm_left_ATL_lin_LASSO_D3_inanimate),10000,2))');
m = mean(group_grOWL_left_inanimate(:,3) - group_lin_LASSO_left_inanimate(:,3));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(3,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_grOWL_inanimate(:,:,1)) - perm_left_ATL_ant_lin_LASSO_D1_inanimate),10000,2))');
m = mean(group_grOWL_left_inanimate(:,4) - group_lin_LASSO_left_inanimate(:,4));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(3,4) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_grOWL_inanimate(:,:,2)) - perm_left_ATL_ant_lin_LASSO_D2_inanimate),10000,2))');
m = mean(group_grOWL_left_inanimate(:,5) - group_lin_LASSO_left_inanimate(:,5));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(3,5) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_ant_grOWL_inanimate(:,:,3)) - perm_left_ATL_ant_lin_LASSO_D3_inanimate),10000,2))');
m = mean(group_grOWL_left_inanimate(:,6) - group_lin_LASSO_left_inanimate(:,6));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(3,6) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_grOWL_inanimate(:,:,1)) - perm_left_ATL_pos_lin_LASSO_D1_inanimate),10000,2))');
m = mean(group_grOWL_left_inanimate(:,7) - group_lin_LASSO_left_inanimate(:,7));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(3,7) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_grOWL_inanimate(:,:,2)) - perm_left_ATL_pos_lin_LASSO_D2_inanimate),10000,2))');
m = mean(group_grOWL_left_inanimate(:,8) - group_lin_LASSO_left_inanimate(:,8));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(3,8) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_left_ATL_pos_grOWL_inanimate(:,:,3)) - perm_left_ATL_pos_lin_LASSO_D3_inanimate),10000,2))');
m = mean(group_grOWL_left_inanimate(:,9) - group_lin_LASSO_left_inanimate(:,9));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(3,9) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_grOWL_inanimate(:,:,1)) - perm_right_ATL_lin_LASSO_D1_inanimate),10000,2))');
m = mean(group_grOWL_right_inanimate(:,1) - group_lin_LASSO_right_inanimate(:,1));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(4,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_grOWL_inanimate(:,:,2)) - perm_right_ATL_lin_LASSO_D2_inanimate),10000,2))');
m = mean(group_grOWL_right_inanimate(:,2) - group_lin_LASSO_right_inanimate(:,2));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(4,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_grOWL_inanimate(:,:,3)) - perm_right_ATL_lin_LASSO_D3_inanimate),10000,2))');
m = mean(group_grOWL_right_inanimate(:,3) - group_lin_LASSO_right_inanimate(:,3));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(4,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_grOWL_inanimate(:,:,1)) - perm_right_ATL_ant_lin_LASSO_D1_inanimate),10000,2))');
m = mean(group_grOWL_right_inanimate(:,4) - group_lin_LASSO_right_inanimate(:,4));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(4,4) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_grOWL_inanimate(:,:,2)) - perm_right_ATL_ant_lin_LASSO_D2_inanimate),10000,2))');
m = mean(group_grOWL_right_inanimate(:,5) - group_lin_LASSO_right_inanimate(:,5));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(4,5) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_ant_grOWL_inanimate(:,:,3)) - perm_right_ATL_ant_lin_LASSO_D3_inanimate),10000,2))');
m = mean(group_grOWL_right_inanimate(:,6) - group_lin_LASSO_right_inanimate(:,6));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(4,6) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_grOWL_inanimate(:,:,1)) - perm_right_ATL_pos_lin_LASSO_D1_inanimate),10000,2))');
m = mean(group_grOWL_right_inanimate(:,7) - group_lin_LASSO_right_inanimate(:,7));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(4,7) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_grOWL_inanimate(:,:,2)) - perm_right_ATL_pos_lin_LASSO_D2_inanimate),10000,2))');
m = mean(group_grOWL_right_inanimate(:,8) - group_lin_LASSO_right_inanimate(:,8));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(4,8) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_right_ATL_pos_grOWL_inanimate(:,:,3)) - perm_right_ATL_pos_lin_LASSO_D3_inanimate),10000,2))');
m = mean(group_grOWL_right_inanimate(:,9) - group_lin_LASSO_right_inanimate(:,9));
b = sum(tmp > m);
inanimate_correlation_difference_p_val(4,9) = (b + 1)/(10000 + 1);


% control the false discovery rate
for i = 1:size(inanimate_correlation_difference_p_val,1)
    p = inanimate_correlation_difference_p_val(i,:);
    inanimate_correlation_difference_p_val(i,:) = mafdr(p,'BHFDR',1);
end

% 2F. Plot results for the ATL (within domain)

% plot animate

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_lin_LASSO_left_animate));
b.FaceColor = 'flat';
b.CData = colours_LASSO_left;
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_lin_LASSO_left_animate)/sqrt(size(group_lin_LASSO_left_animate,1)));
errorbar(1:size(group_lin_LASSO_left_animate,2),mean(group_lin_LASSO_left_animate,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(animate_correlation_p_val,2)
    if (animate_correlation_p_val(1,i) < 0.05 && animate_correlation_difference_p_val(1,i) < 0.05)
       stars{i} = '*†';
    elseif animate_correlation_p_val(1,i) < 0.05
        stars{i} = '*';
    elseif animate_correlation_p_val(1,i) >= 0.05
       stars{i} = '';
    end
end
text((1:9)-0.1,mean(group_lin_LASSO_left_animate)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/lin_LASSO_left_animate.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_lin_LASSO_right_animate));
b.FaceColor = 'flat';
b.CData = colours_LASSO_right;
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_lin_LASSO_right_animate)/sqrt(size(group_lin_LASSO_right_animate,1)));
errorbar(1:size(group_lin_LASSO_right_animate,2),mean(group_lin_LASSO_right_animate,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(animate_correlation_p_val,2)
    if (animate_correlation_p_val(2,i) < 0.05 && animate_correlation_difference_p_val(2,i) < 0.05)
       stars{i} = '*†';
    elseif animate_correlation_p_val(2,i) < 0.05
        stars{i} = '*';
    elseif animate_correlation_p_val(2,i) >= 0.05
       stars{i} = '';
    end
end
text((1:9)-0.1,mean(group_lin_LASSO_right_animate)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/lin_LASSO_right_animate.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_grOWL_left_animate));
b.FaceColor = 'flat';
b.CData = colours_grOWL_left;
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_grOWL_left_animate)/sqrt(size(group_grOWL_left_animate,1)));
errorbar(1:size(group_grOWL_left_animate,2),mean(group_grOWL_left_animate,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(animate_correlation_p_val,2)
    if (animate_correlation_p_val(3,i) < 0.05 && animate_correlation_difference_p_val(3,i) < 0.05)
       stars{i} = '*†';
    elseif animate_correlation_p_val(3,i) < 0.05
        stars{i} = '*';
    elseif animate_correlation_p_val(3,i) >= 0.05
       stars{i} = '';
    end
end
text((1:9)-0.1,mean(group_grOWL_left_animate)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/grOWL_left_animate.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_grOWL_right_animate));
b.FaceColor = 'flat';
b.CData = colours_grOWL_right;
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_grOWL_right_animate)/sqrt(size(group_grOWL_right_animate,1)));
errorbar(1:size(group_grOWL_right_animate,2),mean(group_grOWL_right_animate,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(animate_correlation_p_val,2)
    if (animate_correlation_p_val(4,i) < 0.05 && animate_correlation_difference_p_val(4,i) < 0.05)
       stars{i} = '*†';
    elseif animate_correlation_p_val(4,i) < 0.05
        stars{i} = '*';
    elseif animate_correlation_p_val(4,i) >= 0.05
       stars{i} = '';
    end
end
text((1:9)-0.1,mean(group_grOWL_right_animate)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/grOWL_right_animate.png')
close(fig)

% plot inanimate

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_lin_LASSO_left_inanimate));
b.FaceColor = 'flat';
b.CData = colours_LASSO_left;
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_lin_LASSO_left_inanimate)/sqrt(size(group_lin_LASSO_left_inanimate,1)));
errorbar(1:size(group_lin_LASSO_left_inanimate,2),mean(group_lin_LASSO_left_inanimate,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(inanimate_correlation_p_val,2)
    if (inanimate_correlation_p_val(1,i) < 0.05 && inanimate_correlation_difference_p_val(1,i) < 0.05)
       stars{i} = '*†';
    elseif inanimate_correlation_p_val(1,i) < 0.05
        stars{i} = '*';
    elseif inanimate_correlation_p_val(1,i) >= 0.05
       stars{i} = '';
    end
end
text((1:9)-0.1,mean(group_lin_LASSO_left_inanimate)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/lin_LASSO_left_inanimate.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_lin_LASSO_right_inanimate));
b.FaceColor = 'flat';
b.CData = colours_LASSO_right;
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_lin_LASSO_right_inanimate)/sqrt(size(group_lin_LASSO_right_inanimate,1)));
errorbar(1:size(group_lin_LASSO_right_inanimate,2),mean(group_lin_LASSO_right_inanimate,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(inanimate_correlation_p_val,2)
    if (inanimate_correlation_p_val(2,i) < 0.05 && inanimate_correlation_difference_p_val(2,i) < 0.05)
       stars{i} = '*†';
    elseif inanimate_correlation_p_val(2,i) < 0.05
        stars{i} = '*';
    elseif inanimate_correlation_p_val(2,i) >= 0.05
       stars{i} = '';
    end
end
text((1:9)-0.1,mean(group_lin_LASSO_right_inanimate)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/lin_LASSO_right_inanimate.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_grOWL_left_inanimate));
b.FaceColor = 'flat';
b.CData = colours_grOWL_left;
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_grOWL_left_inanimate)/sqrt(size(group_grOWL_left_inanimate,1)));
errorbar(1:size(group_grOWL_left_inanimate,2),mean(group_grOWL_left_inanimate,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(inanimate_correlation_p_val,2)
    if (inanimate_correlation_p_val(3,i) < 0.05 && inanimate_correlation_difference_p_val(3,i) < 0.05)
       stars{i} = '*†';
    elseif inanimate_correlation_p_val(3,i) < 0.05
        stars{i} = '*';
    elseif inanimate_correlation_p_val(3,i) >= 0.05
       stars{i} = '';
    end
end
text((1:9)-0.1,mean(group_grOWL_left_inanimate)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/grOWL_left_inanimate.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_grOWL_right_inanimate));
b.FaceColor = 'flat';
b.CData = colours_grOWL_right;
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_grOWL_right_inanimate)/sqrt(size(group_grOWL_right_inanimate,1)));
errorbar(1:size(group_grOWL_right_inanimate,2),mean(group_grOWL_right_inanimate,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(inanimate_correlation_p_val,2)
    if (inanimate_correlation_p_val(4,i) < 0.05 && inanimate_correlation_difference_p_val(4,i) < 0.05)
       stars{i} = '*†';
    elseif inanimate_correlation_p_val(4,i) < 0.05
        stars{i} = '*';
    elseif inanimate_correlation_p_val(4,i) >= 0.05
       stars{i} = '';
    end
end
text((1:9)-0.1,mean(group_grOWL_right_inanimate)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/grOWL_right_inanimate.png')
close(fig)

% 3. WHOLE BRAIN (CLASSIFICATION AND CORRELATION)

% 3A. Statistics with Steltzer bootstrapping

% Initialise collated results
group_log_LASSO_wholebrain = zeros(size(subcode,2),1);
group_SOSLASSO_wholebrain = zeros(size(subcode,2),1);
group_lin_LASSO_wholebrain = zeros(size(subcode,2),3);
group_lin_LASSO_wholebrain_animate = zeros(size(subcode,2),3);
group_lin_LASSO_wholebrain_inanimate = zeros(size(subcode,2),3);
group_grOWL_wholebrain = zeros(size(subcode,2),3);
group_grOWL_wholebrain_animate = zeros(size(subcode,2),3);
group_grOWL_wholebrain_inanimate = zeros(size(subcode,2),3);

% collate
group_log_LASSO_wholebrain(:,1) = mean(whole_brain_log_LASSO);
group_SOSLASSO_wholebrain(:,1) = mean(whole_brain_SOSLASSO);
group_lin_LASSO_wholebrain(:,1) = mean(whole_brain_lin_LASSO_D1);
group_lin_LASSO_wholebrain(:,2) = mean(whole_brain_lin_LASSO_D2);
group_lin_LASSO_wholebrain(:,3) = mean(whole_brain_lin_LASSO_D3);
group_lin_LASSO_wholebrain_animate(:,1) = mean(whole_brain_lin_LASSO_D1_animate);
group_lin_LASSO_wholebrain_animate(:,2) = mean(whole_brain_lin_LASSO_D2_animate);
group_lin_LASSO_wholebrain_animate(:,3) = mean(whole_brain_lin_LASSO_D3_animate);
group_lin_LASSO_wholebrain_inanimate(:,1) = mean(whole_brain_lin_LASSO_D1_inanimate);
group_lin_LASSO_wholebrain_inanimate(:,2) = mean(whole_brain_lin_LASSO_D2_inanimate);
group_lin_LASSO_wholebrain_inanimate(:,3) = mean(whole_brain_lin_LASSO_D3_inanimate);
group_grOWL_wholebrain(:,1) = mean(squeeze(whole_brain_grOWL(:,:,1)));
group_grOWL_wholebrain(:,2) = mean(squeeze(whole_brain_grOWL(:,:,2)));
group_grOWL_wholebrain(:,3) = mean(squeeze(whole_brain_grOWL(:,:,3)));
group_grOWL_wholebrain_animate(:,1) = mean(squeeze(whole_brain_grOWL_animate(:,:,1)));
group_grOWL_wholebrain_animate(:,2) = mean(squeeze(whole_brain_grOWL_animate(:,:,2)));
group_grOWL_wholebrain_animate(:,3) = mean(squeeze(whole_brain_grOWL_animate(:,:,3)));
group_grOWL_wholebrain_inanimate(:,1) = mean(squeeze(whole_brain_grOWL_inanimate(:,:,1)));
group_grOWL_wholebrain_inanimate(:,2) = mean(squeeze(whole_brain_grOWL_inanimate(:,:,2)));
group_grOWL_wholebrain_inanimate(:,3) = mean(squeeze(whole_brain_grOWL_inanimate(:,:,3)));

% construct permutation distributions by randomly sampling group averages
% (see Cox et al., 2024, Imaging Neuroscience for methods). Sorting makes
% it easier to caluclate percentile p values
bootstrapped_perm_log_LASSO_wholebrain(:,1) = sort(mean(datasample(perm_whole_brain_log_LASSO,10000,2))');
bootstrapped_perm_SOSLASSO_wholebrain(:,1) = sort(mean(datasample(perm_whole_brain_SOSLASSO,10000,2))');
bootstrapped_perm_lin_LASSO_wholebrain(:,1) = sort(mean(datasample(perm_whole_brain_lin_LASSO_D1,10000,2))');
bootstrapped_perm_lin_LASSO_wholebrain(:,2) = sort(mean(datasample(perm_whole_brain_lin_LASSO_D2,10000,2))');
bootstrapped_perm_lin_LASSO_wholebrain(:,3) = sort(mean(datasample(perm_whole_brain_lin_LASSO_D3,10000,2))');
bootstrapped_perm_lin_LASSO_wholebrain_animate(:,1) = sort(mean(datasample(perm_whole_brain_lin_LASSO_D1_animate,10000,2))');
bootstrapped_perm_lin_LASSO_wholebrain_animate(:,2) = sort(mean(datasample(perm_whole_brain_lin_LASSO_D2_animate,10000,2))');
bootstrapped_perm_lin_LASSO_wholebrain_animate(:,3) = sort(mean(datasample(perm_whole_brain_lin_LASSO_D3_animate,10000,2))');
bootstrapped_perm_lin_LASSO_wholebrain_inanimate(:,1) = sort(mean(datasample(perm_whole_brain_lin_LASSO_D1_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_wholebrain_inanimate(:,2) = sort(mean(datasample(perm_whole_brain_lin_LASSO_D2_inanimate,10000,2))');
bootstrapped_perm_lin_LASSO_wholebrain_inanimate(:,3) = sort(mean(datasample(perm_whole_brain_lin_LASSO_D3_inanimate,10000,2))');
tmp = perm_whole_brain_grOWL(:,:,1);
bootstrapped_perm_grOWL_wholebrain(:,1) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_whole_brain_grOWL(:,:,2);
bootstrapped_perm_grOWL_wholebrain(:,2) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_whole_brain_grOWL(:,:,3);
bootstrapped_perm_grOWL_wholebrain(:,3) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_whole_brain_grOWL_animate(:,:,1);
bootstrapped_perm_grOWL_wholebrain_animate(:,1) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_whole_brain_grOWL_animate(:,:,2);
bootstrapped_perm_grOWL_wholebrain_animate(:,2) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_whole_brain_grOWL_animate(:,:,3);
bootstrapped_perm_grOWL_wholebrain_animate(:,3) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_whole_brain_grOWL_inanimate(:,:,1);
bootstrapped_perm_grOWL_wholebrain_inanimate(:,1) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_whole_brain_grOWL_inanimate(:,:,2);
bootstrapped_perm_grOWL_wholebrain_inanimate(:,2) = sort(mean(datasample(tmp,10000,2))');
tmp = perm_whole_brain_grOWL_inanimate(:,:,3);
bootstrapped_perm_grOWL_wholebrain_inanimate(:,3) = sort(mean(datasample(tmp,10000,2))');

wholebrain_classification_p_val = zeros(2,1);
m = mean(group_log_LASSO_wholebrain(:,1));
b = sum(bootstrapped_perm_log_LASSO_wholebrain(:,1) > m);
wholebrain_classification_p_val(1,1) = (b + 1)/(10000 + 1);
m = mean(group_SOSLASSO_wholebrain(:,1));
b = sum(bootstrapped_perm_SOSLASSO_wholebrain(:,1) > m);
wholebrain_classification_p_val(2,1) = (b + 1)/(10000 + 1);

% control the false discovery rate
for i = 1:size(wholebrain_classification_p_val,1)
    p = wholebrain_classification_p_val(i,:);
    wholebrain_classfication_p_val(i,:) = mafdr(p,'BHFDR',1);
end

% test for differences between methods
wholebrain_classification_difference_p_val = zeros(2,1);
tmp = sort(mean(datasample((perm_whole_brain_log_LASSO - perm_whole_brain_SOSLASSO),10000,2))');
m = mean(group_log_LASSO_wholebrain(:,1) - group_SOSLASSO_wholebrain(:,1));
b = sum(tmp > m);
wholebrain_classification_difference_p_val(1,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_whole_brain_SOSLASSO - perm_whole_brain_log_LASSO),10000,2))');
m = mean(group_SOSLASSO_wholebrain(:,1) - group_log_LASSO_wholebrain(:,1));
b = sum(tmp > m);
wholebrain_classification_difference_p_val(2,1) = (b + 1)/(10000 + 1);
% one test per "row" - no need to control FDR. 

wholebrain_correlation_p_val = zeros(2,3);
m = mean(group_lin_LASSO_wholebrain(:,1));
b = sum(bootstrapped_perm_lin_LASSO_wholebrain(:,1) > m);
wholebrain_correlation_p_val(1,1) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_wholebrain(:,2));
b = sum(bootstrapped_perm_lin_LASSO_wholebrain(:,2) > m);
wholebrain_correlation_p_val(1,2) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_wholebrain(:,3));
b = sum(bootstrapped_perm_lin_LASSO_wholebrain(:,3) > m);
wholebrain_correlation_p_val(1,3) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_wholebrain(:,1));
b = sum(bootstrapped_perm_grOWL_wholebrain(:,1) > m);
wholebrain_correlation_p_val(2,1) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_wholebrain(:,2));
b = sum(bootstrapped_perm_grOWL_wholebrain(:,2) > m);
wholebrain_correlation_p_val(2,2) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_wholebrain(:,3));
b = sum(bootstrapped_perm_grOWL_wholebrain(:,3) > m);
wholebrain_correlation_p_val(2,3) = (b + 1)/(10000 + 1);

% control the false discovery rate
for i = 1:size(wholebrain_correlation_p_val,1)
    p = wholebrain_correlation_p_val(i,:);
    wholebrain_correlation_p_val(i,:) = mafdr(p,'BHFDR',1);
end

% test for differences between methods
wholebrain_correlation_difference_p_val = zeros(2,3);
tmp = sort(mean(datasample(squeeze(perm_whole_brain_lin_LASSO_D1 - perm_whole_brain_grOWL(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_wholebrain(:,1) - group_grOWL_wholebrain(:,1));
b = sum(tmp > m);
wholebrain_correlation_difference_p_val(1,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample(squeeze(perm_whole_brain_lin_LASSO_D2 - perm_whole_brain_grOWL(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_wholebrain(:,2) - group_grOWL_wholebrain(:,2));
b = sum(tmp > m);
wholebrain_correlation_difference_p_val(1,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample(squeeze(perm_whole_brain_lin_LASSO_D3 - perm_whole_brain_grOWL(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_wholebrain(:,3) - group_grOWL_wholebrain(:,3));
b = sum(tmp > m);
wholebrain_correlation_difference_p_val(1,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_whole_brain_grOWL(:,:,1)) - perm_whole_brain_lin_LASSO_D1),10000,2))');
m = mean(group_grOWL_wholebrain(:,1) - group_lin_LASSO_wholebrain(:,1));
b = sum(tmp > m);
wholebrain_correlation_difference_p_val(2,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_whole_brain_grOWL(:,:,2)) - perm_whole_brain_lin_LASSO_D2),10000,2))');
m = mean(group_grOWL_wholebrain(:,2) - group_lin_LASSO_wholebrain(:,2));
b = sum(tmp > m);
wholebrain_correlation_difference_p_val(2,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_whole_brain_grOWL(:,:,3)) - perm_whole_brain_lin_LASSO_D3),10000,2))');
m = mean(group_grOWL_wholebrain(:,3) - group_lin_LASSO_wholebrain(:,3));
b = sum(tmp > m);
wholebrain_correlation_difference_p_val(2,3) = (b + 1)/(10000 + 1);

% control the false discovery rate
for i = 1:size(wholebrain_correlation_difference_p_val,1)
    p = wholebrain_correlation_difference_p_val(i,:);
    wholebrain_correlation_difference_p_val(i,:) = mafdr(p,'BHFDR',1);
end

wholebrain_animate_correlation_p_val = zeros(2,3);
m = mean(group_lin_LASSO_wholebrain_animate(:,1));
b = sum(bootstrapped_perm_lin_LASSO_wholebrain_animate(:,1) > m);
wholebrain_animate_correlation_p_val(1,1) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_wholebrain_animate(:,2));
b = sum(bootstrapped_perm_lin_LASSO_wholebrain_animate(:,2) > m);
wholebrain_animate_correlation_p_val(1,2) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_wholebrain_animate(:,3));
b = sum(bootstrapped_perm_lin_LASSO_wholebrain_animate(:,3) > m);
wholebrain_animate_correlation_p_val(1,3) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_wholebrain_animate(:,1));
b = sum(bootstrapped_perm_grOWL_wholebrain_animate(:,1) > m);
wholebrain_animate_correlation_p_val(2,1) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_wholebrain_animate(:,2));
b = sum(bootstrapped_perm_grOWL_wholebrain_animate(:,2) > m);
wholebrain_animate_correlation_p_val(2,2) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_wholebrain_animate(:,3));
b = sum(bootstrapped_perm_grOWL_wholebrain_animate(:,3) > m);
wholebrain_animate_correlation_p_val(2,3) = (b + 1)/(10000 + 1);

for i = 1:size(wholebrain_animate_correlation_p_val,1)
    p = wholebrain_animate_correlation_p_val(i,:);
    wholebrain_animate_correlation_p_val(i,:) = mafdr(p,'BHFDR',1);
end

% test for differences between methods
wholebrain_animate_correlation_difference_p_val = zeros(2,3);
tmp = sort(mean(datasample(squeeze(perm_whole_brain_lin_LASSO_D1_animate - perm_whole_brain_grOWL_animate(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_wholebrain_animate(:,1) - group_grOWL_wholebrain_animate(:,1));
b = sum(tmp > m);
wholebrain_animate_correlation_difference_p_val(1,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample(squeeze(perm_whole_brain_lin_LASSO_D2_animate - perm_whole_brain_grOWL_animate(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_wholebrain_animate(:,2) - group_grOWL_wholebrain_animate(:,2));
b = sum(tmp > m);
wholebrain_animate_correlation_difference_p_val(1,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample(squeeze(perm_whole_brain_lin_LASSO_D3_animate - perm_whole_brain_grOWL_animate(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_wholebrain_animate(:,3) - group_grOWL_wholebrain_animate(:,3));
b = sum(tmp > m);
wholebrain_animate_correlation_difference_p_val(1,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_whole_brain_grOWL_animate(:,:,1)) - perm_whole_brain_lin_LASSO_D1_animate),10000,2))');
m = mean(group_grOWL_wholebrain_animate(:,1) - group_lin_LASSO_wholebrain_animate(:,1));
b = sum(tmp > m);
wholebrain_animate_correlation_difference_p_val(2,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_whole_brain_grOWL_animate(:,:,2)) - perm_whole_brain_lin_LASSO_D2_animate),10000,2))');
m = mean(group_grOWL_wholebrain_animate(:,2) - group_lin_LASSO_wholebrain_animate(:,2));
b = sum(tmp > m);
wholebrain_animate_correlation_difference_p_val(2,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_whole_brain_grOWL_animate(:,:,3)) - perm_whole_brain_lin_LASSO_D3_animate),10000,2))');
m = mean(group_grOWL_wholebrain_animate(:,3) - group_lin_LASSO_wholebrain_animate(:,3));
b = sum(tmp > m);
wholebrain_animate_correlation_difference_p_val(2,3) = (b + 1)/(10000 + 1);

% control the false discovery rate
for i = 1:size(wholebrain_animate_correlation_difference_p_val,1)
    p = wholebrain_animate_correlation_difference_p_val(i,:);
    wholebrain_animate_correlation_difference_p_val(i,:) = mafdr(p,'BHFDR',1);
end

wholebrain_inanimate_correlation_p_val = zeros(2,3);
m = mean(group_lin_LASSO_wholebrain_inanimate(:,1));
b = sum(bootstrapped_perm_lin_LASSO_wholebrain_inanimate(:,1) > m);
wholebrain_inanimate_correlation_p_val(1,1) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_wholebrain_inanimate(:,2));
b = sum(bootstrapped_perm_lin_LASSO_wholebrain_inanimate(:,2) > m);
wholebrain_inanimate_correlation_p_val(1,2) = (b + 1)/(10000 + 1);
m = mean(group_lin_LASSO_wholebrain_inanimate(:,3));
b = sum(bootstrapped_perm_lin_LASSO_wholebrain_inanimate(:,3) > m);
wholebrain_inanimate_correlation_p_val(1,3) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_wholebrain_inanimate(:,1));
b = sum(bootstrapped_perm_grOWL_wholebrain_inanimate(:,1) > m);
wholebrain_inanimate_correlation_p_val(2,1) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_wholebrain_inanimate(:,2));
b = sum(bootstrapped_perm_grOWL_wholebrain_inanimate(:,2) > m);
wholebrain_inanimate_correlation_p_val(2,2) = (b + 1)/(10000 + 1);
m = mean(group_grOWL_wholebrain_inanimate(:,3));
b = sum(bootstrapped_perm_grOWL_wholebrain_inanimate(:,3) > m);
wholebrain_inanimate_correlation_p_val(2,3) = (b + 1)/(10000 + 1);

for i = 1:size(wholebrain_inanimate_correlation_p_val,1)
    p = wholebrain_inanimate_correlation_p_val(i,:);
    wholebrain_inanimate_correlation_p_val(i,:) = mafdr(p,'BHFDR',1);
end

% test for differences between methods
wholebrain_inanimate_correlation_difference_p_val = zeros(2,3);
tmp = sort(mean(datasample(squeeze(perm_whole_brain_lin_LASSO_D1_inanimate - perm_whole_brain_grOWL_inanimate(:,:,1)),10000,2))');
m = mean(group_lin_LASSO_wholebrain_inanimate(:,1) - group_grOWL_wholebrain_inanimate(:,1));
b = sum(tmp > m);
wholebrain_inanimate_correlation_difference_p_val(1,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample(squeeze(perm_whole_brain_lin_LASSO_D2_inanimate - perm_whole_brain_grOWL_inanimate(:,:,2)),10000,2))');
m = mean(group_lin_LASSO_wholebrain_inanimate(:,2) - group_grOWL_wholebrain_inanimate(:,2));
b = sum(tmp > m);
wholebrain_inanimate_correlation_difference_p_val(1,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample(squeeze(perm_whole_brain_lin_LASSO_D3_inanimate - perm_whole_brain_grOWL_inanimate(:,:,3)),10000,2))');
m = mean(group_lin_LASSO_wholebrain_inanimate(:,3) - group_grOWL_wholebrain_inanimate(:,3));
b = sum(tmp > m);
wholebrain_inanimate_correlation_difference_p_val(1,3) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_whole_brain_grOWL_inanimate(:,:,1)) - perm_whole_brain_lin_LASSO_D1_inanimate),10000,2))');
m = mean(group_grOWL_wholebrain_inanimate(:,1) - group_lin_LASSO_wholebrain_inanimate(:,1));
b = sum(tmp > m);
wholebrain_inanimate_correlation_difference_p_val(2,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_whole_brain_grOWL_inanimate(:,:,2)) - perm_whole_brain_lin_LASSO_D2_inanimate),10000,2))');
m = mean(group_grOWL_wholebrain_inanimate(:,2) - group_lin_LASSO_wholebrain_inanimate(:,2));
b = sum(tmp > m);
wholebrain_inanimate_correlation_difference_p_val(2,2) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((squeeze(perm_whole_brain_grOWL_inanimate(:,:,3)) - perm_whole_brain_lin_LASSO_D3_inanimate),10000,2))');
m = mean(group_grOWL_wholebrain_inanimate(:,3) - group_lin_LASSO_wholebrain_inanimate(:,3));
b = sum(tmp > m);
wholebrain_inanimate_correlation_difference_p_val(2,3) = (b + 1)/(10000 + 1);

% control the false discovery rate
for i = 1:size(wholebrain_inanimate_correlation_difference_p_val,1)
    p = wholebrain_inanimate_correlation_difference_p_val(i,:);
    wholebrain_inanimate_correlation_difference_p_val(i,:) = mafdr(p,'BHFDR',1);
end

% 3B. Plot

% plot classification

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_log_LASSO_wholebrain));
b.FaceColor = 'flat';
b.CData = colours_LASSO_wholebrain(1,:);
set(gca,'xticklabel',{})
ylabel('Accuracy')
set(gca,'FontSize',20)
ylim([0.4,1])
yline(0.5,'--')
hold on
confint = 1.96*(std(group_log_LASSO_wholebrain)/sqrt(size(group_log_LASSO_wholebrain,1)));
errorbar(1:size(group_log_LASSO_wholebrain,2),mean(group_log_LASSO_wholebrain,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(wholebrain_classification_p_val,2)
    if (wholebrain_classification_p_val(1,i) < 0.05 && wholebrain_classification_difference_p_val(1,i) < 0.05)
       stars{i} = '*†';
    elseif wholebrain_classification_p_val(1,i) < 0.05
        stars{i} = '*';
    elseif wholebrain_classification_p_val(1,i) >= 0.05
       stars{i} = '';
    end
end
% fix for stars going off the page
text(0.9,0.975,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/log_LASSO_wholebrain.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_SOSLASSO_wholebrain));
b.FaceColor = 'flat';
b.CData = colours_SOSLASSO_wholebrain(1,:);
set(gca,'xticklabel',{})
ylabel('Accuracy')
set(gca,'FontSize',20)
ylim([0.4,1])
yline(0.5,'--')
hold on
confint = 1.96*(std(group_SOSLASSO_wholebrain)/sqrt(size(group_SOSLASSO_wholebrain,1)));
errorbar(1:size(group_SOSLASSO_wholebrain,2),mean(group_SOSLASSO_wholebrain,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(wholebrain_classification_p_val,2)
    if (wholebrain_classification_p_val(2,i) < 0.05 && wholebrain_classification_difference_p_val(2,i) < 0.05)
       stars{i} = '*†';
    elseif wholebrain_classification_p_val(2,i) < 0.05
        stars{i} = '*';
    elseif wholebrain_classification_p_val(2,i) >= 0.05
       stars{i} = '';
    end
end
% fix for stars going off the page
text(0.9,0.975,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/SOSLASSO_wholebrain.png')
close(fig)

% plot correlation

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_lin_LASSO_wholebrain));
b.FaceColor = 'flat';
b.CData = [colours_LASSO_wholebrain(1,:);colours_LASSO_wholebrain(1,:);colours_LASSO_wholebrain(1,:)];
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_lin_LASSO_wholebrain)/sqrt(size(group_lin_LASSO_wholebrain,1)));
errorbar(1:size(group_lin_LASSO_wholebrain,2),mean(group_lin_LASSO_wholebrain,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(wholebrain_correlation_p_val,2)
    if (wholebrain_correlation_p_val(1,i) < 0.05 && wholebrain_correlation_difference_p_val(1,i) < 0.05)
       stars{i} = '*†';
    elseif wholebrain_correlation_p_val(1,i) < 0.05
        stars{i} = '*';
    elseif wholebrain_correlation_p_val(1,i) >= 0.05
       stars{i} = '';
    end
end
text((1:3)-0.1,mean(group_lin_LASSO_wholebrain)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/lin_LASSO_wholebrain.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_grOWL_wholebrain));
b.FaceColor = 'flat';
b.CData = [colours_grOWL_wholebrain(1,:);colours_grOWL_wholebrain(1,:);colours_grOWL_wholebrain(1,:)];
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_grOWL_wholebrain)/sqrt(size(group_grOWL_wholebrain,1)));
errorbar(1:size(group_grOWL_wholebrain,2),mean(group_grOWL_wholebrain,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(wholebrain_correlation_p_val,2)
    if (wholebrain_correlation_p_val(2,i) < 0.05 && wholebrain_correlation_difference_p_val(2,i) < 0.05)
       stars{i} = '*†';
    elseif wholebrain_correlation_p_val(2,i) < 0.05
        stars{i} = '*';
    elseif wholebrain_correlation_p_val(2,i) >= 0.05
       stars{i} = '';
    end
end
text((1:3)-0.1,mean(group_grOWL_wholebrain)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/grOWL_wholebrain.png')
close(fig)

% plot animate
fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_lin_LASSO_wholebrain_animate));
b.FaceColor = 'flat';
b.CData = [colours_LASSO_wholebrain(1,:);colours_LASSO_wholebrain(1,:);colours_LASSO_wholebrain(1,:)];
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_lin_LASSO_wholebrain_animate)/sqrt(size(group_lin_LASSO_wholebrain_animate,1)));
errorbar(1:size(group_lin_LASSO_wholebrain_animate,2),mean(group_lin_LASSO_wholebrain_animate,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(wholebrain_animate_correlation_p_val,2)
    if (wholebrain_animate_correlation_p_val(1,i) < 0.05 && wholebrain_animate_correlation_difference_p_val(1,i) < 0.05)
       stars{i} = '*†';
    elseif wholebrain_animate_correlation_p_val(1,i) < 0.05
        stars{i} = '*';
    elseif wholebrain_animate_correlation_p_val(1,i) >= 0.05
       stars{i} = '';
    end
end
text((1:3)-0.1,mean(group_lin_LASSO_wholebrain_animate)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/lin_LASSO_wholebrain_animate.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_grOWL_wholebrain_animate));
b.FaceColor = 'flat';
b.CData = [colours_grOWL_wholebrain(1,:);colours_grOWL_wholebrain(1,:);colours_grOWL_wholebrain(1,:)];
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_grOWL_wholebrain_animate)/sqrt(size(group_grOWL_wholebrain_animate,1)));
errorbar(1:size(group_grOWL_wholebrain_animate,2),mean(group_grOWL_wholebrain_animate,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(wholebrain_animate_correlation_p_val,2)
    if (wholebrain_animate_correlation_p_val(2,i) < 0.05 && wholebrain_animate_correlation_difference_p_val(2,i) < 0.05)
       stars{i} = '*†';
    elseif wholebrain_animate_correlation_p_val(2,i) < 0.05
        stars{i} = '*';
    elseif wholebrain_animate_correlation_p_val(2,i) >= 0.05
       stars{i} = '';
    end
end
text((1:3)-0.1,mean(group_grOWL_wholebrain_animate)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/grOWL_wholebrain_animate.png')
close(fig)

% plot inanimate

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_lin_LASSO_wholebrain_inanimate));
b.FaceColor = 'flat';
b.CData = [colours_LASSO_wholebrain(1,:);colours_LASSO_wholebrain(1,:);colours_LASSO_wholebrain(1,:)];
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_lin_LASSO_wholebrain_inanimate)/sqrt(size(group_lin_LASSO_wholebrain_inanimate,1)));
errorbar(1:size(group_lin_LASSO_wholebrain_inanimate,2),mean(group_lin_LASSO_wholebrain_inanimate,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(wholebrain_inanimate_correlation_p_val,2)
    if (wholebrain_inanimate_correlation_p_val(1,i) < 0.05 && wholebrain_inanimate_correlation_difference_p_val(1,i) < 0.05)
       stars{i} = '*†';
    elseif wholebrain_inanimate_correlation_p_val(1,i) < 0.05
        stars{i} = '*';
    elseif wholebrain_inanimate_correlation_p_val(1,i) >= 0.05
       stars{i} = '';
    end
end
text((1:3)-0.1,mean(group_lin_LASSO_wholebrain_inanimate)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/lin_LASSO_wholebrain_inanimate.png')
close(fig)

fig = figure;
fig.Position = [610 240 900 725];
b = bar(mean(group_grOWL_wholebrain_inanimate));
b.FaceColor = 'flat';
b.CData = [colours_grOWL_wholebrain(1,:);colours_grOWL_wholebrain(1,:);colours_grOWL_wholebrain(1,:)];
set(gca,'xticklabel',{'D1','D2','D3','D1','D2','D3','D1','D2','D3'})
ylabel('Correlation')
set(gca,'FontSize',20)
ylim([-0.1,1])
hold on
confint = 1.96*(std(group_grOWL_wholebrain_inanimate)/sqrt(size(group_grOWL_wholebrain_inanimate,1)));
errorbar(1:size(group_grOWL_wholebrain_inanimate,2),mean(group_grOWL_wholebrain_inanimate,1),confint,'LineStyle','none','Color','black');
stars = {};
for i =1:size(wholebrain_inanimate_correlation_p_val,2)
    if (wholebrain_inanimate_correlation_p_val(2,i) < 0.05 && wholebrain_inanimate_correlation_difference_p_val(2,i) < 0.05)
       stars{i} = '*†';
    elseif wholebrain_inanimate_correlation_p_val(2,i) < 0.05
        stars{i} = '*';
    elseif wholebrain_inanimate_correlation_p_val(2,i) >= 0.05
       stars{i} = '';
    end
end
text((1:3)-0.1,mean(group_grOWL_wholebrain_inanimate)+confint+0.05,stars,'FontSize',20)
saveas(fig,'/group/mlr-lab/Saskia/7T_decoding/bar_charts/grOWL_wholebrain_inanimate.png')
close(fig)


%% test for dynamism (animacy only)

dynamism_p_val = zeros(4,1);
tmp = sort(mean(datasample((perm_left_ATL_pos_log_LASSO - perm_left_ATL_ant_log_LASSO),10000,2))');
m = mean(group_log_LASSO_left(:,3) - group_log_LASSO_left(:,2));
b = sum(tmp > m);
dynamism_p_val(1,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_right_ATL_pos_log_LASSO - perm_right_ATL_ant_log_LASSO),10000,2))');
m = mean(group_log_LASSO_right(:,3) - group_log_LASSO_right(:,2));
b = sum(tmp > m);
dynamism_p_val(2,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_left_ATL_pos_SOSLASSO - perm_left_ATL_ant_SOSLASSO),10000,2))');
m = mean(group_SOSLASSO_left(:,3) - group_SOSLASSO_left(:,2));
b = sum(tmp > m);
dynamism_p_val(3,1) = (b + 1)/(10000 + 1);
tmp = sort(mean(datasample((perm_right_ATL_pos_SOSLASSO - perm_right_ATL_ant_SOSLASSO),10000,2))');
m = mean(group_SOSLASSO_right(:,3) - group_SOSLASSO_right(:,2));
b = sum(tmp > m);
dynamism_p_val(4,1) = (b + 1)/(10000 + 1);

