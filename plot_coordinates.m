% makes graphs of coordinates (target and predicted)

addpath('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts');

% target

% load labels
load('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/stimulimaster.mat');
labels = cellstr(stimulimaster);
% load similarity matrix
load('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/dilkina_norms.mat');
% decompose similarity matrix
U = embed_similarity_matrix(dilkina_norms,3);

% set up colours
% categorise (1 - land mammals, 2 - birds, 3 - insects, 4 - amphibians, 5 -
% reptiles, 6 - marine animals, 7 - musical instruments, 8 - vehicles, 9 -
% clothing, 10 - other inanimate natural, 11 - other inanimate manmade, 12
% - buildings)
cidx = [1, 2, 6, 5, 3, 5, 6, 3, 3, 1, 3, 1, 1, 3, 1, 1, 1, 2, 2, 1, 6, 3, 1, 4, 1, 1, 1, 1, 1, 1, 6, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 6, 1, 3, 5, 3, 1, 2, 1, 1, 7, 10, 11, 11, 8, 10, 11, 8, 11, 11, 11, 7, 11, 11, 11, 12, 11, 11, 9, 9, 7, 11, 7, 9, 7, 11, 11, 7, 8, 11, 11, 11, 9, 8, 11, 7, 11, 9, 11, 10, 11, 11, 11, 8, 7,	11,	7, 11, 11, 12];
% assign each item a colour
c = orderedcolors('gem12');
colours = zeros(length(stimulimaster),3);
for i = 1:length(stimulimaster)
    colours(i,:) = c(cidx(i),:);
end
% plot - D1 v. D2
figure;
subplot(2,2,1);
p = scatter(U(:,1),U(:,2),'filled');
p.CData = colours;
xlabel('Dimension 1')
ylabel('Dimension 2')
hold on
% generate slight offset for labels
dx = 0.001;
dy = 0.001;
text(U(:,1)+dx,U(:,2)+dy,labels);
% D1 v. D3
subplot(2,2,2);
p = scatter(U(:,1),U(:,3),'filled');
p.CData = colours;
xlabel('Dimension 1')
ylabel('Dimension 3')
hold on
text(U(:,1)+dx,U(:,3)+dy,labels);

% plot 3d, just for fun
figure;
p = scatter3(U(:,1),U(:,2),U(:,3),'filled');
p.CData = colours;
xlabel('Dimension 1')
ylabel('Dimension 2')
zlabel('Dimension 3')
hold on
dz = 0.001;
text(U(:,1)+dx,U(:,2)+dy,U(:,3)+dz,labels);
