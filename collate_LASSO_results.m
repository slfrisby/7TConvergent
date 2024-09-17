% by Saskia. Loads, collates and plots LASSO results generated with
% fit_models. 

root = ['/imaging/projects/cbu/wbic-p00591-DAISY/main/derivatives'];
cd([root]);

subcode = {'sub-001','sub-002','sub-003','sub-004','sub-007','sub-009','sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019','sub-020','sub-021','sub-022','sub-023','sub-024','sub-026','sub-028','sub-029','sub-030','sub-031','sub-032'};


% initialise variables
% across-run decoding is 40-fold
acrossRun_whole_brain_t2star = zeros(40,size(subcode,2));
acrossRun_left_ATL_t2star = zeros(40,size(subcode,2));
acrossRun_right_ATL_t2star = zeros(40,size(subcode,2));
acrossRun_left_ATL_ant_t2star = zeros(40,size(subcode,2));
acrossRun_left_ATL_pos_t2star = zeros(40,size(subcode,2));
acrossRun_right_ATL_ant_t2star = zeros(40,size(subcode,2));
acrossRun_right_ATL_pos_t2star = zeros(40,size(subcode,2));
acrossRun_whole_brain_tedana = zeros(40,size(subcode,2));
acrossRun_left_ATL_tedana = zeros(40,size(subcode,2));
acrossRun_right_ATL_tedana = zeros(40,size(subcode,2));
acrossRun_left_ATL_ant_tedana = zeros(40,size(subcode,2));
acrossRun_left_ATL_pos_tedana = zeros(40,size(subcode,2));
acrossRun_right_ATL_ant_tedana = zeros(40,size(subcode,2));
acrossRun_right_ATL_pos_tedana = zeros(40,size(subcode,2));
% decoding on averaged data is 10-fold 
averaged_whole_brain_t2star = zeros(10,size(subcode,2));
averaged_left_ATL_t2star = zeros(10,size(subcode,2));
averaged_right_ATL_t2star = zeros(10,size(subcode,2));
averaged_left_ATL_ant_t2star = zeros(10,size(subcode,2));
averaged_left_ATL_pos_t2star = zeros(10,size(subcode,2));
averaged_right_ATL_ant_t2star = zeros(10,size(subcode,2));
averaged_right_ATL_pos_t2star = zeros(10,size(subcode,2));
averaged_whole_brain_tedana = zeros(10,size(subcode,2));
averaged_left_ATL_tedana = zeros(10,size(subcode,2));
averaged_right_ATL_tedana = zeros(10,size(subcode,2));
averaged_left_ATL_ant_tedana = zeros(10,size(subcode,2));
averaged_left_ATL_pos_tedana = zeros(10,size(subcode,2));
averaged_right_ATL_ant_tedana = zeros(10,size(subcode,2));
averaged_right_ATL_pos_tedana = zeros(10,size(subcode,2));

% for every participant
for s = 1:size(subcode,2)
        
    % load t2star results
    result = load([root,'/LASSO/t2star/',subcode{s},'/results.mat']);

    % add results to collated matrix
    acrossRun_whole_brain_t2star(:,s) = result.acrossRun_whole_brain;
    acrossRun_left_ATL_t2star(:,s) = result.acrossRun_left_ATL;
    acrossRun_right_ATL_t2star(:,s) = result.acrossRun_right_ATL;
    acrossRun_left_ATL_ant_t2star(:,s) = result.acrossRun_left_ATL_ant;
    acrossRun_left_ATL_pos_t2star(:,s) = result.acrossRun_left_ATL_pos;
    acrossRun_right_ATL_ant_t2star(:,s) = result.acrossRun_right_ATL_ant;
    acrossRun_right_ATL_pos_t2star(:,s) = result.acrossRun_right_ATL_pos;
    averaged_whole_brain_t2star(:,s) = result.averaged_whole_brain;
    averaged_left_ATL_t2star(:,s) = result.averaged_left_ATL;
    averaged_right_ATL_t2star(:,s) = result.averaged_right_ATL;
    averaged_left_ATL_ant_t2star(:,s) = result.averaged_left_ATL_ant;
    averaged_left_ATL_pos_t2star(:,s) = result.averaged_left_ATL_pos;
    averaged_right_ATL_ant_t2star(:,s) = result.averaged_right_ATL_ant;
    averaged_right_ATL_pos_t2star(:,s) = result.averaged_right_ATL_pos;

    % load tedana results
    result = load([root,'/LASSO/tedana/',subcode{s},'/results.mat']);

    % add results to collated matrix
    acrossRun_whole_brain_tedana(:,s) = result.acrossRun_whole_brain;
    acrossRun_left_ATL_tedana(:,s) = result.acrossRun_left_ATL;
    acrossRun_right_ATL_tedana(:,s) = result.acrossRun_right_ATL;
    acrossRun_left_ATL_ant_tedana(:,s) = result.acrossRun_left_ATL_ant;
    acrossRun_left_ATL_pos_tedana(:,s) = result.acrossRun_left_ATL_pos;
    acrossRun_right_ATL_ant_tedana(:,s) = result.acrossRun_right_ATL_ant;
    acrossRun_right_ATL_pos_tedana(:,s) = result.acrossRun_right_ATL_pos;
    averaged_whole_brain_tedana(:,s) = result.averaged_whole_brain;
    averaged_left_ATL_tedana(:,s) = result.averaged_left_ATL;
    averaged_right_ATL_tedana(:,s) = result.averaged_right_ATL;
    averaged_left_ATL_ant_tedana(:,s) = result.averaged_left_ATL_ant;
    averaged_left_ATL_pos_tedana(:,s) = result.averaged_left_ATL_pos;
    averaged_right_ATL_ant_tedana(:,s) = result.averaged_right_ATL_ant;
    averaged_right_ATL_pos_tedana(:,s) = result.averaged_right_ATL_pos;

end

% save everything in a single .mat file
save([root,'/LASSO/all_results.mat'],'acrossRun_whole_brain_t2star','acrossRun_left_ATL_t2star','acrossRun_right_ATL_t2star','acrossRun_left_ATL_ant_t2star','acrossRun_left_ATL_pos_t2star','acrossRun_right_ATL_ant_t2star','acrossRun_right_ATL_ant_t2star','acrossRun_right_ATL_pos_t2star','averaged_whole_brain_t2star','averaged_left_ATL_t2star','averaged_right_ATL_t2star','averaged_left_ATL_ant_t2star','averaged_left_ATL_pos_t2star','averaged_right_ATL_ant_t2star','averaged_right_ATL_pos_t2star','acrossRun_whole_brain_tedana','acrossRun_left_ATL_tedana','acrossRun_right_ATL_tedana','acrossRun_left_ATL_ant_tedana','acrossRun_left_ATL_pos_tedana','acrossRun_right_ATL_ant_tedana','acrossRun_right_ATL_pos_tedana','averaged_whole_brain_tedana','averaged_left_ATL_tedana','averaged_right_ATL_tedana','averaged_left_ATL_ant_tedana','averaged_left_ATL_pos_tedana','averaged_right_ATL_ant_tedana','averaged_right_ATL_pos_tedana');

% plot! First figure shows performance for every participant
figure;
% across run - whole brain - t2star
subplot(6,2,1)
bar(subcode,mean(acrossRun_whole_brain_t2star,1),'FaceColor','#77AC30')
ylim([0.4,1])
yline(0.5,'--')
title('Across run - whole brain - t2star')
hold on
confint = 1.96*(std(acrossRun_whole_brain_t2star)/sqrt(size(acrossRun_whole_brain_t2star,1)));
errorbar(1:size(subcode,2),mean(acrossRun_whole_brain_t2star,1),confint,'LineStyle','none','Color','black');
% across run - whole brain - tedana
subplot(6,2,3)
bar(subcode,mean(acrossRun_whole_brain_tedana,1),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([0.4,1])
yline(0.5,'--')
title('Across run - whole brain - tedana')
hold on
confint = 1.96*(std(acrossRun_whole_brain_tedana)/sqrt(size(acrossRun_whole_brain_tedana,1)));
errorbar(1:size(subcode,2),mean(acrossRun_whole_brain_tedana,1),confint,'LineStyle','none','Color','black');
% across run - left ATL - t2star
subplot(6,2,5)
bar(subcode,mean(acrossRun_left_ATL_t2star,1),'FaceColor','#A2142F')
ylim([0.4,1])
yline(0.5,'--')
title('Across run - left ATL - t2star')
hold on
confint = 1.96*(std(acrossRun_left_ATL_t2star)/sqrt(size(acrossRun_left_ATL_t2star,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_t2star,1),confint,'LineStyle','none','Color','black');
% across run - left ATL - tedana
subplot(6,2,7)
bar(subcode,mean(acrossRun_left_ATL_tedana,1),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([0.4,1])
yline(0.5,'--')
title('Across run - left ATL - tedana')
hold on
confint = 1.96*(std(acrossRun_left_ATL_tedana)/sqrt(size(acrossRun_left_ATL_tedana,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_tedana,1),confint,'LineStyle','none','Color','black');
% across run - right ATL - t2star
subplot(6,2,9)
bar(subcode,mean(acrossRun_right_ATL_t2star,1),'FaceColor','#0072BD')
ylim([0.4,1])
yline(0.5,'--')
title('Across run - right ATL - t2star')
hold on
confint = 1.96*(std(acrossRun_right_ATL_t2star)/sqrt(size(acrossRun_right_ATL_t2star,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_t2star,1),confint,'LineStyle','none','Color','black');
% across run - right ATL - tedana
subplot(6,2,11)
bar(subcode,mean(acrossRun_right_ATL_tedana,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([0.4,1])
yline(0.5,'--')
title('Across run - right ATL - tedana')
hold on
confint = 1.96*(std(acrossRun_right_ATL_tedana)/sqrt(size(acrossRun_right_ATL_tedana,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_tedana,1),confint,'LineStyle','none','Color','black');
% averaged - whole brain - t2star
subplot(6,2,2)
bar(subcode,mean(averaged_whole_brain_t2star,1),'FaceColor','#77AC30')
ylim([0.4,1])
yline(0.5,'--')
title('Averaged - whole brain - t2star')
hold on
confint = 1.96*(std(averaged_whole_brain_t2star)/sqrt(size(averaged_whole_brain_t2star,1)));
errorbar(1:size(subcode,2),mean(averaged_whole_brain_t2star,1),confint,'LineStyle','none','Color','black');
% averaged - whole brain - tedana
subplot(6,2,4)
bar(subcode,mean(averaged_whole_brain_tedana,1),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([0.4,1])
yline(0.5,'--')
title('Averaged - whole brain - tedana')
hold on
confint = 1.96*(std(averaged_whole_brain_tedana)/sqrt(size(averaged_whole_brain_tedana,1)));
errorbar(1:size(subcode,2),mean(averaged_whole_brain_tedana,1),confint,'LineStyle','none','Color','black');
% averaged - left ATL - t2star
subplot(6,2,6)
bar(subcode,mean(averaged_left_ATL_t2star,1),'FaceColor','#A2142F')
ylim([0.4,1])
yline(0.5,'--')
title('Averaged - left ATL - t2star')
hold on
confint = 1.96*(std(averaged_left_ATL_t2star)/sqrt(size(averaged_left_ATL_t2star,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_t2star,1),confint,'LineStyle','none','Color','black');
% averaged - left ATL - tedana
subplot(6,2,8)
bar(subcode,mean(averaged_left_ATL_tedana,1),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([0.4,1])
yline(0.5,'--')
title('Averaged - left ATL - tedana')
hold on
confint = 1.96*(std(averaged_left_ATL_tedana)/sqrt(size(averaged_left_ATL_tedana,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_tedana,1),confint,'LineStyle','none','Color','black');
% averaged - right ATL - t2star
subplot(6,2,10)
bar(subcode,mean(averaged_right_ATL_t2star,1),'FaceColor','#0072BD')
ylim([0.4,1])
yline(0.5,'--')
title('Averaged - right ATL - t2star')
hold on
confint = 1.96*(std(averaged_right_ATL_t2star)/sqrt(size(averaged_right_ATL_t2star,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_t2star,1),confint,'LineStyle','none','Color','black');
% averaged - right ATL - tedana
subplot(6,2,12)
bar(subcode,mean(averaged_right_ATL_tedana,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([0.4,1])
yline(0.5,'--')
title('Averaged - right ATL - tedana')
hold on
confint = 1.96*(std(averaged_right_ATL_tedana)/sqrt(size(averaged_right_ATL_tedana,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_tedana,1),confint,'LineStyle','none','Color','black');
sgtitle('LASSO results')

% second figure shows group averages
group = zeros(size(subcode,2),12);
group(:,1) = mean(acrossRun_whole_brain_t2star)';
group(:,2) = mean(acrossRun_whole_brain_tedana)';
group(:,3) = mean(acrossRun_left_ATL_t2star)';
group(:,4) = mean(acrossRun_left_ATL_tedana)';
group(:,5) = mean(acrossRun_right_ATL_t2star)';
group(:,6) = mean(acrossRun_right_ATL_tedana)';
group(:,7) = mean(averaged_whole_brain_t2star)';
group(:,8) = mean(averaged_whole_brain_tedana)';
group(:,9) = mean(averaged_left_ATL_t2star)';
group(:,10) = mean(averaged_left_ATL_tedana)';
group(:,11) = mean(averaged_right_ATL_t2star)';
group(:,12) = mean(averaged_right_ATL_tedana)';
colnames = {'AR-WB-t2star','AR-WB-tedana','AR-lATL-t2star','AR-lATL-tedana','AR-rATL-t2star','AR-rATL-tedana','av-WB-t2star','av-WB-tedana','av-lATL-t2star','av-lATL-tedana','av-rATL-t2star','av-rATL-tedana'};

figure;
b = bar(colnames,mean(group),'FaceColor','flat');
% normalised RGB triplets
colours = [119,172,48;119,172,48;162,20,47;162,20,47;0,114,189;0,114,189;119,172,48;119,172,48;162,20,47;162,20,47;0,114,189;0,114,189]/255;
b.CData = colours;
ylim([0.4,1])
yline(0.5,'--')
title('Group results')
hold on
confint = 1.96*(std(group)/sqrt(size(group,1)));
errorbar(1:size(group,2),mean(group,1),confint,'LineStyle','none','Color','black');

% third figure shows results for anterior and posterior ROIs for every
% participant
figure;
% across run - left ROI (anterior) - t2star
subplot(8,2,1)
bar(subcode,mean(acrossRun_left_ATL_ant_t2star,1),'FaceColor','#D95319')
ylim([0.4,1])
yline(0.5,'--')
title('Across run - left ROI (anterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_left_ATL_ant_t2star)/sqrt(size(acrossRun_left_ATL_ant_t2star,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_ant_t2star,1),confint,'LineStyle','none','Color','black');
% across run - left ROI (anterior) - tedana
subplot(8,2,3)
bar(subcode,mean(acrossRun_left_ATL_ant_tedana,1),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([0.4,1])
yline(0.5,'--')
title('Across run - left ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_left_ATL_ant_tedana)/sqrt(size(acrossRun_left_ATL_ant_tedana,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_ant_tedana,1),confint,'LineStyle','none','Color','black');
% across run - left ROI (posterior) - t2star
subplot(8,2,5)
bar(subcode,mean(acrossRun_left_ATL_pos_t2star,1),'FaceColor','#EDB120')
ylim([0.4,1])
yline(0.5,'--')
title('Across run - left ROI (posterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_left_ATL_pos_t2star)/sqrt(size(acrossRun_left_ATL_pos_t2star,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_pos_t2star,1),confint,'LineStyle','none','Color','black');
% across run - left ROI (posterior) - tedana
subplot(8,2,7)
bar(subcode,mean(acrossRun_left_ATL_pos_tedana,1),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([0.4,1])
yline(0.5,'--')
title('Across run - left ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_left_ATL_pos_tedana)/sqrt(size(acrossRun_left_ATL_pos_tedana,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_pos_tedana,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (anterior) - t2star
subplot(8,2,9)
bar(subcode,mean(acrossRun_right_ATL_ant_t2star,1),'FaceColor','#0072BD')
ylim([0.4,1])
yline(0.5,'--')
title('Across run - right ROI (anterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_right_ATL_ant_t2star)/sqrt(size(acrossRun_right_ATL_ant_t2star,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_ant_t2star,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (anterior) - tedana
subplot(8,2,11)
bar(subcode,mean(acrossRun_right_ATL_ant_tedana,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([0.4,1])
yline(0.5,'--')
title('Across run - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_right_ATL_ant_tedana)/sqrt(size(acrossRun_right_ATL_ant_tedana,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_ant_tedana,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (posterior) - t2star
subplot(8,2,13)
bar(subcode,mean(acrossRun_right_ATL_pos_t2star,1),'FaceColor','#4DBEEE')
ylim([0.4,1])
yline(0.5,'--')
title('Across run - right ROI (posterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_right_ATL_pos_t2star)/sqrt(size(acrossRun_right_ATL_pos_t2star,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_pos_t2star,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (posterior) - tedana
subplot(8,2,15)
bar(subcode,mean(acrossRun_right_ATL_pos_tedana,1),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([0.4,1])
yline(0.5,'--')
title('Across run - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_right_ATL_pos_tedana)/sqrt(size(acrossRun_right_ATL_pos_tedana,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_pos_tedana,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (anterior) - t2star
subplot(8,2,2)
bar(subcode,mean(averaged_left_ATL_ant_t2star,1),'FaceColor','#D95319')
ylim([0.4,1])
yline(0.5,'--')
title('Averaged - left ROI (anterior) - t2star')
hold on
confint = 1.96*(std(averaged_left_ATL_ant_t2star)/sqrt(size(averaged_left_ATL_ant_t2star,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_ant_t2star,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (anterior) - tedana
subplot(8,2,4)
bar(subcode,mean(averaged_left_ATL_ant_tedana,1),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([0.4,1])
yline(0.5,'--')
title('Averaged - left ROI (anterior) - tedana')
hold on
confint = 1.96*(std(averaged_left_ATL_ant_tedana)/sqrt(size(averaged_left_ATL_ant_tedana,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_ant_tedana,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (posterior) - t2star
subplot(8,2,6)
bar(subcode,mean(averaged_left_ATL_pos_t2star,1),'FaceColor','#EDB120')
ylim([0.4,1])
yline(0.5,'--')
title('Averaged - left ROI (posterior) - t2star')
hold on
confint = 1.96*(std(averaged_left_ATL_pos_t2star)/sqrt(size(averaged_left_ATL_pos_t2star,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_pos_t2star,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (posterior) - tedana
subplot(8,2,8)
bar(subcode,mean(averaged_left_ATL_pos_tedana,1),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([0.4,1])
yline(0.5,'--')
title('Averaged - left ROI (anterior) - tedana')
hold on
confint = 1.96*(std(averaged_left_ATL_pos_tedana)/sqrt(size(averaged_left_ATL_pos_tedana,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_pos_tedana,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (anterior) - t2star
subplot(8,2,10)
bar(subcode,mean(averaged_right_ATL_ant_t2star,1),'FaceColor','#0072BD')
ylim([0.4,1])
yline(0.5,'--')
title('Averaged - right ROI (anterior) - t2star')
hold on
confint = 1.96*(std(averaged_right_ATL_ant_t2star)/sqrt(size(averaged_right_ATL_ant_t2star,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_ant_t2star,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (anterior) - tedana
subplot(8,2,12)
bar(subcode,mean(averaged_right_ATL_ant_tedana,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([0.4,1])
yline(0.5,'--')
title('Averaged - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(averaged_right_ATL_ant_tedana)/sqrt(size(averaged_right_ATL_ant_tedana,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_ant_tedana,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (posterior) - t2star
subplot(8,2,14)
bar(subcode,mean(averaged_right_ATL_pos_t2star,1),'FaceColor','#4DBEEE')
ylim([0.4,1])
yline(0.5,'--')
title('Averaged - right ROI (posterior) - t2star')
hold on
confint = 1.96*(std(averaged_right_ATL_pos_t2star)/sqrt(size(averaged_right_ATL_pos_t2star,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_pos_t2star,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (posterior) - tedana
subplot(8,2,16)
bar(subcode,mean(averaged_right_ATL_pos_tedana,1),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([0.4,1])
yline(0.5,'--')
title('Averaged - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(averaged_right_ATL_pos_tedana)/sqrt(size(averaged_right_ATL_pos_tedana,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_pos_tedana,1),confint,'LineStyle','none','Color','black');
sgtitle('LASSO results - anterior and posterior halves of ROI')

% fourth figure shows group averages
group = zeros(size(subcode,2),16);
group(:,1) = mean(acrossRun_left_ATL_ant_t2star)';
group(:,2) = mean(acrossRun_left_ATL_ant_tedana)';
group(:,3) = mean(acrossRun_left_ATL_pos_t2star)';
group(:,4) = mean(acrossRun_left_ATL_pos_tedana)';
group(:,5) = mean(acrossRun_right_ATL_ant_t2star)';
group(:,6) = mean(acrossRun_right_ATL_ant_tedana)';
group(:,7) = mean(acrossRun_right_ATL_pos_t2star)';
group(:,8) = mean(acrossRun_right_ATL_pos_tedana)';
group(:,9) = mean(averaged_left_ATL_ant_t2star)';
group(:,10) = mean(averaged_left_ATL_ant_tedana)';
group(:,11) = mean(averaged_left_ATL_pos_t2star)';
group(:,12) = mean(averaged_left_ATL_pos_tedana)';
group(:,13) = mean(averaged_right_ATL_ant_t2star)';
group(:,14) = mean(averaged_right_ATL_ant_tedana)';
group(:,15) = mean(averaged_right_ATL_pos_t2star)';
group(:,16) = mean(averaged_right_ATL_pos_tedana)';
colnames = {'AR-lATL-ant-t2star','AR-lATL-ant-tedana','AR-lATL-pos-t2star','AR-lATL-pos-tedana','AR-rATL-ant-t2star','AR-rATL-ant-tedana','AR-rATL-pos-t2star','AR-rATL-pos-tedana','av-lATL-ant-t2star','av-lATL-ant-tedana','av-lATL-pos-t2star','av-lATL-pos-tedana','av-rATL-ant-t2star','av-rATL-ant-tedana','av-rATL-pos-t2star','av-rATL-pos-tedana'};

figure;
b = bar(colnames,mean(group),'FaceColor','flat');
% normalised RGB triplets
colours = [0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0,0.4470,0.7410;0,0.4470,0.7410;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0,0.4470,0.7410;0,0.4470,0.7410;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330];
b.CData = colours;
ylim([0.4,1])
yline(0.5,'--')
title('Group results')
hold on
confint = 1.96*(std(group)/sqrt(size(group,1)));
errorbar(1:size(group,2),mean(group,1),confint,'LineStyle','none','Color','black');
