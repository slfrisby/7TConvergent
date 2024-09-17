% by Saskia. Loads, collates and plots linear results generated with
% fit__linear_models. 

root = ['/imaging/projects/cbu/wbic-p00591-DAISY/main/derivatives'];
cd([root]);

subcode = {'sub-001','sub-002','sub-003','sub-004','sub-007','sub-009','sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019','sub-020','sub-021','sub-022','sub-023','sub-024','sub-026','sub-028','sub-029','sub-030','sub-031','sub-032'};

%% collate results

% initialise variables
% overall hold-out accuracy first. Across-run decoding is 40-fold
acrossRun_whole_brain_t2star_D1 = zeros(40,size(subcode,2));
acrossRun_whole_brain_t2star_D2 = zeros(40,size(subcode,2));
acrossRun_whole_brain_t2star_D3 = zeros(40,size(subcode,2));
acrossRun_left_ATL_t2star_D1 = zeros(40,size(subcode,2));
acrossRun_left_ATL_t2star_D2 = zeros(40,size(subcode,2));
acrossRun_left_ATL_t2star_D3 = zeros(40,size(subcode,2));
acrossRun_right_ATL_t2star_D1 = zeros(40,size(subcode,2));
acrossRun_right_ATL_t2star_D2 = zeros(40,size(subcode,2));
acrossRun_right_ATL_t2star_D3 = zeros(40,size(subcode,2));
acrossRun_left_ATL_ant_t2star_D1 = zeros(40,size(subcode,2));
acrossRun_left_ATL_ant_t2star_D2 = zeros(40,size(subcode,2));
acrossRun_left_ATL_ant_t2star_D3 = zeros(40,size(subcode,2));
acrossRun_left_ATL_pos_t2star_D1 = zeros(40,size(subcode,2));
acrossRun_left_ATL_pos_t2star_D2 = zeros(40,size(subcode,2));
acrossRun_left_ATL_pos_t2star_D3 = zeros(40,size(subcode,2));
acrossRun_right_ATL_ant_t2star_D1 = zeros(40,size(subcode,2));
acrossRun_right_ATL_ant_t2star_D2 = zeros(40,size(subcode,2));
acrossRun_right_ATL_ant_t2star_D3 = zeros(40,size(subcode,2));
acrossRun_right_ATL_pos_t2star_D1 = zeros(40,size(subcode,2));
acrossRun_right_ATL_pos_t2star_D2 = zeros(40,size(subcode,2));
acrossRun_right_ATL_pos_t2star_D3 = zeros(40,size(subcode,2));
acrossRun_whole_brain_tedana_D1 = zeros(40,size(subcode,2));
acrossRun_whole_brain_tedana_D2 = zeros(40,size(subcode,2));
acrossRun_whole_brain_tedana_D3 = zeros(40,size(subcode,2));
acrossRun_left_ATL_tedana_D1 = zeros(40,size(subcode,2));
acrossRun_left_ATL_tedana_D2 = zeros(40,size(subcode,2));
acrossRun_left_ATL_tedana_D3 = zeros(40,size(subcode,2));
acrossRun_right_ATL_tedana_D1 = zeros(40,size(subcode,2));
acrossRun_right_ATL_tedana_D2 = zeros(40,size(subcode,2));
acrossRun_right_ATL_tedana_D3 = zeros(40,size(subcode,2));
acrossRun_left_ATL_ant_tedana_D1 = zeros(40,size(subcode,2));
acrossRun_left_ATL_ant_tedana_D2 = zeros(40,size(subcode,2));
acrossRun_left_ATL_ant_tedana_D3 = zeros(40,size(subcode,2));
acrossRun_left_ATL_pos_tedana_D1 = zeros(40,size(subcode,2));
acrossRun_left_ATL_pos_tedana_D2 = zeros(40,size(subcode,2));
acrossRun_left_ATL_pos_tedana_D3 = zeros(40,size(subcode,2));
acrossRun_right_ATL_ant_tedana_D1 = zeros(40,size(subcode,2));
acrossRun_right_ATL_ant_tedana_D2 = zeros(40,size(subcode,2));
acrossRun_right_ATL_ant_tedana_D3 = zeros(40,size(subcode,2));
acrossRun_right_ATL_pos_tedana_D1 = zeros(40,size(subcode,2));
acrossRun_right_ATL_pos_tedana_D2 = zeros(40,size(subcode,2));
acrossRun_right_ATL_pos_tedana_D3 = zeros(40,size(subcode,2));
% decoding on averaged data is 10-fold
averaged_whole_brain_t2star_D1 = zeros(10,size(subcode,2));
averaged_whole_brain_t2star_D2 = zeros(10,size(subcode,2));
averaged_whole_brain_t2star_D3 = zeros(10,size(subcode,2));
averaged_left_ATL_t2star_D1 = zeros(10,size(subcode,2));
averaged_left_ATL_t2star_D2 = zeros(10,size(subcode,2));
averaged_left_ATL_t2star_D3 = zeros(10,size(subcode,2));
averaged_right_ATL_t2star_D1 = zeros(10,size(subcode,2));
averaged_right_ATL_t2star_D2 = zeros(10,size(subcode,2));
averaged_right_ATL_t2star_D3 = zeros(10,size(subcode,2));
averaged_left_ATL_ant_t2star_D1 = zeros(10,size(subcode,2));
averaged_left_ATL_ant_t2star_D2 = zeros(10,size(subcode,2));
averaged_left_ATL_ant_t2star_D3 = zeros(10,size(subcode,2));
averaged_left_ATL_pos_t2star_D1 = zeros(10,size(subcode,2));
averaged_left_ATL_pos_t2star_D2 = zeros(10,size(subcode,2));
averaged_left_ATL_pos_t2star_D3 = zeros(10,size(subcode,2));
averaged_right_ATL_ant_t2star_D1 = zeros(10,size(subcode,2));
averaged_right_ATL_ant_t2star_D2 = zeros(10,size(subcode,2));
averaged_right_ATL_ant_t2star_D3 = zeros(10,size(subcode,2));
averaged_right_ATL_pos_t2star_D1 = zeros(10,size(subcode,2));
averaged_right_ATL_pos_t2star_D2 = zeros(10,size(subcode,2));
averaged_right_ATL_pos_t2star_D3 = zeros(10,size(subcode,2));
averaged_whole_brain_tedana_D1 = zeros(10,size(subcode,2));
averaged_whole_brain_tedana_D2 = zeros(10,size(subcode,2));
averaged_whole_brain_tedana_D3 = zeros(10,size(subcode,2));
averaged_left_ATL_tedana_D1 = zeros(10,size(subcode,2));
averaged_left_ATL_tedana_D2 = zeros(10,size(subcode,2));
averaged_left_ATL_tedana_D3 = zeros(10,size(subcode,2));
averaged_right_ATL_tedana_D1 = zeros(10,size(subcode,2));
averaged_right_ATL_tedana_D2 = zeros(10,size(subcode,2));
averaged_right_ATL_tedana_D3 = zeros(10,size(subcode,2));
averaged_left_ATL_ant_tedana_D1 = zeros(10,size(subcode,2));
averaged_left_ATL_ant_tedana_D2 = zeros(10,size(subcode,2));
averaged_left_ATL_ant_tedana_D3 = zeros(10,size(subcode,2));
averaged_left_ATL_pos_tedana_D1 = zeros(10,size(subcode,2));
averaged_left_ATL_pos_tedana_D2 = zeros(10,size(subcode,2));
averaged_left_ATL_pos_tedana_D3 = zeros(10,size(subcode,2));
averaged_right_ATL_ant_tedana_D1 = zeros(10,size(subcode,2));
averaged_right_ATL_ant_tedana_D2 = zeros(10,size(subcode,2));
averaged_right_ATL_ant_tedana_D3 = zeros(10,size(subcode,2));
averaged_right_ATL_pos_tedana_D1 = zeros(10,size(subcode,2));
averaged_right_ATL_pos_tedana_D2 = zeros(10,size(subcode,2));
averaged_right_ATL_pos_tedana_D3 = zeros(10,size(subcode,2));
% then predicted coordinates. Across-run will estimate predicted
% coordinates 4 times
coords_acrossRun_whole_brain_t2star_D1 = zeros(400,size(subcode,2));
coords_acrossRun_whole_brain_t2star_D2 = zeros(400,size(subcode,2));
coords_acrossRun_whole_brain_t2star_D3 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_t2star_D1 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_t2star_D2 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_t2star_D3 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_t2star_D1 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_t2star_D2 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_t2star_D3 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_ant_t2star_D1 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_ant_t2star_D2 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_ant_t2star_D3 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_pos_t2star_D1 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_pos_t2star_D2 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_pos_t2star_D3 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_ant_t2star_D1 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_ant_t2star_D2 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_ant_t2star_D3 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_pos_t2star_D1 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_pos_t2star_D2 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_pos_t2star_D3 = zeros(400,size(subcode,2));
coords_acrossRun_whole_brain_tedana_D1 = zeros(400,size(subcode,2));
coords_acrossRun_whole_brain_tedana_D2 = zeros(400,size(subcode,2));
coords_acrossRun_whole_brain_tedana_D3 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_tedana_D1 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_tedana_D2 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_tedana_D3 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_tedana_D1 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_tedana_D2 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_tedana_D3 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_ant_tedana_D1 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_ant_tedana_D2 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_ant_tedana_D3 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_pos_tedana_D1 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_pos_tedana_D2 = zeros(400,size(subcode,2));
coords_acrossRun_left_ATL_pos_tedana_D3 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_ant_tedana_D1 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_ant_tedana_D2 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_ant_tedana_D3 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_pos_tedana_D1 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_pos_tedana_D2 = zeros(400,size(subcode,2));
coords_acrossRun_right_ATL_pos_tedana_D3 = zeros(400,size(subcode,2));
% Averaged will estimate predicted coordinates once
coords_averaged_whole_brain_t2star_D1 = zeros(100,size(subcode,2));
coords_averaged_whole_brain_t2star_D2 = zeros(100,size(subcode,2));
coords_averaged_whole_brain_t2star_D3 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_t2star_D1 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_t2star_D2 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_t2star_D3 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_t2star_D1 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_t2star_D2 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_t2star_D3 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_ant_t2star_D1 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_ant_t2star_D2 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_ant_t2star_D3 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_pos_t2star_D1 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_pos_t2star_D2 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_pos_t2star_D3 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_ant_t2star_D1 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_ant_t2star_D2 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_ant_t2star_D3 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_pos_t2star_D1 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_pos_t2star_D2 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_pos_t2star_D3 = zeros(100,size(subcode,2));
coords_averaged_whole_brain_tedana_D1 = zeros(100,size(subcode,2));
coords_averaged_whole_brain_tedana_D2 = zeros(100,size(subcode,2));
coords_averaged_whole_brain_tedana_D3 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_tedana_D1 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_tedana_D2 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_tedana_D3 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_tedana_D1 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_tedana_D2 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_tedana_D3 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_ant_tedana_D1 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_ant_tedana_D2 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_ant_tedana_D3 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_pos_tedana_D1 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_pos_tedana_D2 = zeros(100,size(subcode,2));
coords_averaged_left_ATL_pos_tedana_D3 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_ant_tedana_D1 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_ant_tedana_D2 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_ant_tedana_D3 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_pos_tedana_D1 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_pos_tedana_D2 = zeros(100,size(subcode,2));
coords_averaged_right_ATL_pos_tedana_D3 = zeros(100,size(subcode,2));

% for every participant
for s = 1:size(subcode,2)
        
    % load t2star results
    result = load([root,'/linear/t2star/',subcode{s},'/results.mat']);

    % add results to collated matrix
    acrossRun_whole_brain_t2star_D1(:,s) = result.acrossRun_whole_brain_D1.output;
    acrossRun_whole_brain_t2star_D2(:,s) = result.acrossRun_whole_brain_D2.output;
    acrossRun_whole_brain_t2star_D3(:,s) = result.acrossRun_whole_brain_D3.output;
    acrossRun_left_ATL_t2star_D1(:,s) = result.acrossRun_left_ATL_D1.output;
    acrossRun_left_ATL_t2star_D2(:,s) = result.acrossRun_left_ATL_D2.output;
    acrossRun_left_ATL_t2star_D3(:,s) = result.acrossRun_left_ATL_D3.output;
    acrossRun_right_ATL_t2star_D1(:,s) = result.acrossRun_right_ATL_D1.output;
    acrossRun_right_ATL_t2star_D2(:,s) = result.acrossRun_right_ATL_D2.output;
    acrossRun_right_ATL_t2star_D3(:,s) = result.acrossRun_right_ATL_D3.output;
    acrossRun_left_ATL_ant_t2star_D1(:,s) = result.acrossRun_left_ATL_ant_D1.output;
    acrossRun_left_ATL_ant_t2star_D2(:,s) = result.acrossRun_left_ATL_ant_D2.output;
    acrossRun_left_ATL_ant_t2star_D3(:,s) = result.acrossRun_left_ATL_ant_D3.output;
    acrossRun_left_ATL_pos_t2star_D1(:,s) = result.acrossRun_left_ATL_pos_D1.output;
    acrossRun_left_ATL_pos_t2star_D2(:,s) = result.acrossRun_left_ATL_pos_D2.output;
    acrossRun_left_ATL_pos_t2star_D3(:,s) = result.acrossRun_left_ATL_pos_D3.output;
    acrossRun_right_ATL_ant_t2star_D1(:,s) = result.acrossRun_right_ATL_ant_D1.output;
    acrossRun_right_ATL_ant_t2star_D2(:,s) = result.acrossRun_right_ATL_ant_D2.output;
    acrossRun_right_ATL_ant_t2star_D3(:,s) = result.acrossRun_right_ATL_ant_D3.output;
    acrossRun_right_ATL_pos_t2star_D1(:,s) = result.acrossRun_right_ATL_pos_D1.output;
    acrossRun_right_ATL_pos_t2star_D2(:,s) = result.acrossRun_right_ATL_pos_D2.output;
    acrossRun_right_ATL_pos_t2star_D3(:,s) = result.acrossRun_right_ATL_pos_D3.output;
    averaged_whole_brain_t2star_D1(:,s) = result.averaged_whole_brain_D1.output;
    averaged_whole_brain_t2star_D2(:,s) = result.averaged_whole_brain_D2.output;
    averaged_whole_brain_t2star_D3(:,s) = result.averaged_whole_brain_D3.output;
    averaged_left_ATL_t2star_D1(:,s) = result.averaged_left_ATL_D1.output;
    averaged_left_ATL_t2star_D2(:,s) = result.averaged_left_ATL_D2.output;
    averaged_left_ATL_t2star_D3(:,s) = result.averaged_left_ATL_D3.output;
    averaged_right_ATL_t2star_D1(:,s) = result.averaged_right_ATL_D1.output;
    averaged_right_ATL_t2star_D2(:,s) = result.averaged_right_ATL_D2.output;
    averaged_right_ATL_t2star_D3(:,s) = result.averaged_right_ATL_D3.output;
    averaged_left_ATL_ant_t2star_D1(:,s) = result.averaged_left_ATL_ant_D1.output;
    averaged_left_ATL_ant_t2star_D2(:,s) = result.averaged_left_ATL_ant_D2.output;
    averaged_left_ATL_ant_t2star_D3(:,s) = result.averaged_left_ATL_ant_D3.output;
    averaged_left_ATL_pos_t2star_D1(:,s) = result.averaged_left_ATL_pos_D1.output;
    averaged_left_ATL_pos_t2star_D2(:,s) = result.averaged_left_ATL_pos_D2.output;
    averaged_left_ATL_pos_t2star_D3(:,s) = result.averaged_left_ATL_pos_D3.output;
    averaged_right_ATL_ant_t2star_D1(:,s) = result.averaged_right_ATL_ant_D1.output;
    averaged_right_ATL_ant_t2star_D2(:,s) = result.averaged_right_ATL_ant_D2.output;
    averaged_right_ATL_ant_t2star_D3(:,s) = result.averaged_right_ATL_ant_D3.output;
    averaged_right_ATL_pos_t2star_D1(:,s) = result.averaged_right_ATL_pos_D1.output;
    averaged_right_ATL_pos_t2star_D2(:,s) = result.averaged_right_ATL_pos_D2.output;
    averaged_right_ATL_pos_t2star_D3(:,s) = result.averaged_right_ATL_pos_D3.output;
    % also collate coordinates
    coords_acrossRun_whole_brain_t2star_D1(:,s) = result.acrossRun_whole_brain_D1.predictedcoords;
    coords_acrossRun_whole_brain_t2star_D2(:,s) = result.acrossRun_whole_brain_D2.predictedcoords;
    coords_acrossRun_whole_brain_t2star_D3(:,s) = result.acrossRun_whole_brain_D3.predictedcoords;
    coords_acrossRun_left_ATL_t2star_D1(:,s) = result.acrossRun_left_ATL_D1.predictedcoords;
    coords_acrossRun_left_ATL_t2star_D2(:,s) = result.acrossRun_left_ATL_D2.predictedcoords;
    coords_acrossRun_left_ATL_t2star_D3(:,s) = result.acrossRun_left_ATL_D3.predictedcoords;
    coords_acrossRun_right_ATL_t2star_D1(:,s) = result.acrossRun_right_ATL_D1.predictedcoords;
    coords_acrossRun_right_ATL_t2star_D2(:,s) = result.acrossRun_right_ATL_D2.predictedcoords;
    coords_acrossRun_right_ATL_t2star_D3(:,s) = result.acrossRun_right_ATL_D3.predictedcoords;
    coords_acrossRun_left_ATL_ant_t2star_D1(:,s) = result.acrossRun_left_ATL_ant_D1.predictedcoords;
    coords_acrossRun_left_ATL_ant_t2star_D2(:,s) = result.acrossRun_left_ATL_ant_D2.predictedcoords;
    coords_acrossRun_left_ATL_ant_t2star_D3(:,s) = result.acrossRun_left_ATL_ant_D3.predictedcoords;
    coords_acrossRun_left_ATL_pos_t2star_D1(:,s) = result.acrossRun_left_ATL_pos_D1.predictedcoords;
    coords_acrossRun_left_ATL_pos_t2star_D2(:,s) = result.acrossRun_left_ATL_pos_D2.predictedcoords;
    coords_acrossRun_left_ATL_pos_t2star_D3(:,s) = result.acrossRun_left_ATL_pos_D3.predictedcoords;
    coords_acrossRun_right_ATL_ant_t2star_D1(:,s) = result.acrossRun_right_ATL_ant_D1.predictedcoords;
    coords_acrossRun_right_ATL_ant_t2star_D2(:,s) = result.acrossRun_right_ATL_ant_D2.predictedcoords;
    coords_acrossRun_right_ATL_ant_t2star_D3(:,s) = result.acrossRun_right_ATL_ant_D3.predictedcoords;
    coords_acrossRun_right_ATL_pos_t2star_D1(:,s) = result.acrossRun_right_ATL_pos_D1.predictedcoords;
    coords_acrossRun_right_ATL_pos_t2star_D2(:,s) = result.acrossRun_right_ATL_pos_D2.predictedcoords;
    coords_acrossRun_right_ATL_pos_t2star_D3(:,s) = result.acrossRun_right_ATL_pos_D3.predictedcoords;
    coords_averaged_whole_brain_t2star_D1(:,s) = result.averaged_whole_brain_D1.predictedcoords;
    coords_averaged_whole_brain_t2star_D2(:,s) = result.averaged_whole_brain_D2.predictedcoords;
    coords_averaged_whole_brain_t2star_D3(:,s) = result.averaged_whole_brain_D3.predictedcoords;
    coords_averaged_left_ATL_t2star_D1(:,s) = result.averaged_left_ATL_D1.predictedcoords;
    coords_averaged_left_ATL_t2star_D2(:,s) = result.averaged_left_ATL_D2.predictedcoords;
    coords_averaged_left_ATL_t2star_D3(:,s) = result.averaged_left_ATL_D3.predictedcoords;
    coords_averaged_right_ATL_t2star_D1(:,s) = result.averaged_right_ATL_D1.predictedcoords;
    coords_averaged_right_ATL_t2star_D2(:,s) = result.averaged_right_ATL_D2.predictedcoords;
    coords_averaged_right_ATL_t2star_D3(:,s) = result.averaged_right_ATL_D3.predictedcoords;
    coords_averaged_left_ATL_ant_t2star_D1(:,s) = result.averaged_left_ATL_ant_D1.predictedcoords;
    coords_averaged_left_ATL_ant_t2star_D2(:,s) = result.averaged_left_ATL_ant_D2.predictedcoords;
    coords_averaged_left_ATL_ant_t2star_D3(:,s) = result.averaged_left_ATL_ant_D3.predictedcoords;
    coords_averaged_left_ATL_pos_t2star_D1(:,s) = result.averaged_left_ATL_pos_D1.predictedcoords;
    coords_averaged_left_ATL_pos_t2star_D2(:,s) = result.averaged_left_ATL_pos_D2.predictedcoords;
    coords_averaged_left_ATL_pos_t2star_D3(:,s) = result.averaged_left_ATL_pos_D3.predictedcoords;
    coords_averaged_right_ATL_ant_t2star_D1(:,s) = result.averaged_right_ATL_ant_D1.predictedcoords;
    coords_averaged_right_ATL_ant_t2star_D2(:,s) = result.averaged_right_ATL_ant_D2.predictedcoords;
    coords_averaged_right_ATL_ant_t2star_D3(:,s) = result.averaged_right_ATL_ant_D3.predictedcoords;
    coords_averaged_right_ATL_pos_t2star_D1(:,s) = result.averaged_right_ATL_pos_D1.predictedcoords;
    coords_averaged_right_ATL_pos_t2star_D2(:,s) = result.averaged_right_ATL_pos_D2.predictedcoords;
    coords_averaged_right_ATL_pos_t2star_D3(:,s) = result.averaged_right_ATL_pos_D3.predictedcoords;

    % load tedana results
    result = load([root,'/linear/tedana/',subcode{s},'/results.mat']);

    % add results to collated matrix
    acrossRun_whole_brain_tedana_D1(:,s) = result.acrossRun_whole_brain_D1.output;
    acrossRun_whole_brain_tedana_D2(:,s) = result.acrossRun_whole_brain_D2.output;
    acrossRun_whole_brain_tedana_D3(:,s) = result.acrossRun_whole_brain_D3.output;
    acrossRun_left_ATL_tedana_D1(:,s) = result.acrossRun_left_ATL_D1.output;
    acrossRun_left_ATL_tedana_D2(:,s) = result.acrossRun_left_ATL_D2.output;
    acrossRun_left_ATL_tedana_D3(:,s) = result.acrossRun_left_ATL_D3.output;
    acrossRun_right_ATL_tedana_D1(:,s) = result.acrossRun_right_ATL_D1.output;
    acrossRun_right_ATL_tedana_D2(:,s) = result.acrossRun_right_ATL_D2.output;
    acrossRun_right_ATL_tedana_D3(:,s) = result.acrossRun_right_ATL_D3.output;
    acrossRun_left_ATL_ant_tedana_D1(:,s) = result.acrossRun_left_ATL_ant_D1.output;
    acrossRun_left_ATL_ant_tedana_D2(:,s) = result.acrossRun_left_ATL_ant_D2.output;
    acrossRun_left_ATL_ant_tedana_D3(:,s) = result.acrossRun_left_ATL_ant_D3.output;
    acrossRun_left_ATL_pos_tedana_D1(:,s) = result.acrossRun_left_ATL_pos_D1.output;
    acrossRun_left_ATL_pos_tedana_D2(:,s) = result.acrossRun_left_ATL_pos_D2.output;
    acrossRun_left_ATL_pos_tedana_D3(:,s) = result.acrossRun_left_ATL_pos_D3.output;
    acrossRun_right_ATL_ant_tedana_D1(:,s) = result.acrossRun_right_ATL_ant_D1.output;
    acrossRun_right_ATL_ant_tedana_D2(:,s) = result.acrossRun_right_ATL_ant_D2.output;
    acrossRun_right_ATL_ant_tedana_D3(:,s) = result.acrossRun_right_ATL_ant_D3.output;
    acrossRun_right_ATL_pos_tedana_D1(:,s) = result.acrossRun_right_ATL_pos_D1.output;
    acrossRun_right_ATL_pos_tedana_D2(:,s) = result.acrossRun_right_ATL_pos_D2.output;
    acrossRun_right_ATL_pos_tedana_D3(:,s) = result.acrossRun_right_ATL_pos_D3.output;
    averaged_whole_brain_tedana_D1(:,s) = result.averaged_whole_brain_D1.output;
    averaged_whole_brain_tedana_D2(:,s) = result.averaged_whole_brain_D2.output;
    averaged_whole_brain_tedana_D3(:,s) = result.averaged_whole_brain_D3.output;
    averaged_left_ATL_tedana_D1(:,s) = result.averaged_left_ATL_D1.output;
    averaged_left_ATL_tedana_D2(:,s) = result.averaged_left_ATL_D2.output;
    averaged_left_ATL_tedana_D3(:,s) = result.averaged_left_ATL_D3.output;
    averaged_right_ATL_tedana_D1(:,s) = result.averaged_right_ATL_D1.output;
    averaged_right_ATL_tedana_D2(:,s) = result.averaged_right_ATL_D2.output;
    averaged_right_ATL_tedana_D3(:,s) = result.averaged_right_ATL_D3.output;
    averaged_left_ATL_ant_tedana_D1(:,s) = result.averaged_left_ATL_ant_D1.output;
    averaged_left_ATL_ant_tedana_D2(:,s) = result.averaged_left_ATL_ant_D2.output;
    averaged_left_ATL_ant_tedana_D3(:,s) = result.averaged_left_ATL_ant_D3.output;
    averaged_left_ATL_pos_tedana_D1(:,s) = result.averaged_left_ATL_pos_D1.output;
    averaged_left_ATL_pos_tedana_D2(:,s) = result.averaged_left_ATL_pos_D2.output;
    averaged_left_ATL_pos_tedana_D3(:,s) = result.averaged_left_ATL_pos_D3.output;
    averaged_right_ATL_ant_tedana_D1(:,s) = result.averaged_right_ATL_ant_D1.output;
    averaged_right_ATL_ant_tedana_D2(:,s) = result.averaged_right_ATL_ant_D2.output;
    averaged_right_ATL_ant_tedana_D3(:,s) = result.averaged_right_ATL_ant_D3.output;
    averaged_right_ATL_pos_tedana_D1(:,s) = result.averaged_right_ATL_pos_D1.output;
    averaged_right_ATL_pos_tedana_D2(:,s) = result.averaged_right_ATL_pos_D2.output;
    averaged_right_ATL_pos_tedana_D3(:,s) = result.averaged_right_ATL_pos_D3.output;
    % also collate coordinates
    coords_acrossRun_whole_brain_tedana_D1(:,s) = result.acrossRun_whole_brain_D1.predictedcoords;
    coords_acrossRun_whole_brain_tedana_D2(:,s) = result.acrossRun_whole_brain_D2.predictedcoords;
    coords_acrossRun_whole_brain_tedana_D3(:,s) = result.acrossRun_whole_brain_D3.predictedcoords;
    coords_acrossRun_left_ATL_tedana_D1(:,s) = result.acrossRun_left_ATL_D1.predictedcoords;
    coords_acrossRun_left_ATL_tedana_D2(:,s) = result.acrossRun_left_ATL_D2.predictedcoords;
    coords_acrossRun_left_ATL_tedana_D3(:,s) = result.acrossRun_left_ATL_D3.predictedcoords;
    coords_acrossRun_right_ATL_tedana_D1(:,s) = result.acrossRun_right_ATL_D1.predictedcoords;
    coords_acrossRun_right_ATL_tedana_D2(:,s) = result.acrossRun_right_ATL_D2.predictedcoords;
    coords_acrossRun_right_ATL_tedana_D3(:,s) = result.acrossRun_right_ATL_D3.predictedcoords;
    coords_acrossRun_left_ATL_ant_tedana_D1(:,s) = result.acrossRun_left_ATL_ant_D1.predictedcoords;
    coords_acrossRun_left_ATL_ant_tedana_D2(:,s) = result.acrossRun_left_ATL_ant_D2.predictedcoords;
    coords_acrossRun_left_ATL_ant_tedana_D3(:,s) = result.acrossRun_left_ATL_ant_D3.predictedcoords;
    coords_acrossRun_left_ATL_pos_tedana_D1(:,s) = result.acrossRun_left_ATL_pos_D1.predictedcoords;
    coords_acrossRun_left_ATL_pos_tedana_D2(:,s) = result.acrossRun_left_ATL_pos_D2.predictedcoords;
    coords_acrossRun_left_ATL_pos_tedana_D3(:,s) = result.acrossRun_left_ATL_pos_D3.predictedcoords;
    coords_acrossRun_right_ATL_ant_tedana_D1(:,s) = result.acrossRun_right_ATL_ant_D1.predictedcoords;
    coords_acrossRun_right_ATL_ant_tedana_D2(:,s) = result.acrossRun_right_ATL_ant_D2.predictedcoords;
    coords_acrossRun_right_ATL_ant_tedana_D3(:,s) = result.acrossRun_right_ATL_ant_D3.predictedcoords;
    coords_acrossRun_right_ATL_pos_tedana_D1(:,s) = result.acrossRun_right_ATL_pos_D1.predictedcoords;
    coords_acrossRun_right_ATL_pos_tedana_D2(:,s) = result.acrossRun_right_ATL_pos_D2.predictedcoords;
    coords_acrossRun_right_ATL_pos_tedana_D3(:,s) = result.acrossRun_right_ATL_pos_D3.predictedcoords;
    coords_averaged_whole_brain_tedana_D1(:,s) = result.averaged_whole_brain_D1.predictedcoords;
    coords_averaged_whole_brain_tedana_D2(:,s) = result.averaged_whole_brain_D2.predictedcoords;
    coords_averaged_whole_brain_tedana_D3(:,s) = result.averaged_whole_brain_D3.predictedcoords;
    coords_averaged_left_ATL_tedana_D1(:,s) = result.averaged_left_ATL_D1.predictedcoords;
    coords_averaged_left_ATL_tedana_D2(:,s) = result.averaged_left_ATL_D2.predictedcoords;
    coords_averaged_left_ATL_tedana_D3(:,s) = result.averaged_left_ATL_D3.predictedcoords;
    coords_averaged_right_ATL_tedana_D1(:,s) = result.averaged_right_ATL_D1.predictedcoords;
    coords_averaged_right_ATL_tedana_D2(:,s) = result.averaged_right_ATL_D2.predictedcoords;
    coords_averaged_right_ATL_tedana_D3(:,s) = result.averaged_right_ATL_D3.predictedcoords;
    coords_averaged_left_ATL_ant_tedana_D1(:,s) = result.averaged_left_ATL_ant_D1.predictedcoords;
    coords_averaged_left_ATL_ant_tedana_D2(:,s) = result.averaged_left_ATL_ant_D2.predictedcoords;
    coords_averaged_left_ATL_ant_tedana_D3(:,s) = result.averaged_left_ATL_ant_D3.predictedcoords;
    coords_averaged_left_ATL_pos_tedana_D1(:,s) = result.averaged_left_ATL_pos_D1.predictedcoords;
    coords_averaged_left_ATL_pos_tedana_D2(:,s) = result.averaged_left_ATL_pos_D2.predictedcoords;
    coords_averaged_left_ATL_pos_tedana_D3(:,s) = result.averaged_left_ATL_pos_D3.predictedcoords;
    coords_averaged_right_ATL_ant_tedana_D1(:,s) = result.averaged_right_ATL_ant_D1.predictedcoords;
    coords_averaged_right_ATL_ant_tedana_D2(:,s) = result.averaged_right_ATL_ant_D2.predictedcoords;
    coords_averaged_right_ATL_ant_tedana_D3(:,s) = result.averaged_right_ATL_ant_D3.predictedcoords;
    coords_averaged_right_ATL_pos_tedana_D1(:,s) = result.averaged_right_ATL_pos_D1.predictedcoords;
    coords_averaged_right_ATL_pos_tedana_D2(:,s) = result.averaged_right_ATL_pos_D2.predictedcoords;
    coords_averaged_right_ATL_pos_tedana_D3(:,s) = result.averaged_right_ATL_pos_D3.predictedcoords;

    clear result

end

% save everything in a single .mat file
save([root,'/linear/all_results.mat'])

%% overall accuracy across folds

% plot! First figure shows dimension 1 for every participant
figure;
% across run - whole brain - t2star
subplot(6,2,1)
bar(subcode,mean(acrossRun_whole_brain_t2star_D1,1),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Across run - whole brain - t2star')
hold on
confint = 1.96*(std(acrossRun_whole_brain_t2star_D1)/sqrt(size(acrossRun_whole_brain_t2star_D1,1)));
errorbar(1:size(subcode,2),mean(acrossRun_whole_brain_t2star_D1,1),confint,'LineStyle','none','Color','black');
% across run - whole brain - tedana
subplot(6,2,3)
bar(subcode,mean(acrossRun_whole_brain_tedana_D1,1),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - whole brain - tedana')
hold on
confint = 1.96*(std(acrossRun_whole_brain_tedana_D1)/sqrt(size(acrossRun_whole_brain_tedana_D1,1)));
errorbar(1:size(subcode,2),mean(acrossRun_whole_brain_tedana_D1,1),confint,'LineStyle','none','Color','black');
% across run - left ATL - t2star
subplot(6,2,5)
bar(subcode,mean(acrossRun_left_ATL_t2star_D1,1),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Across run - left ATL - t2star')
hold on
confint = 1.96*(std(acrossRun_left_ATL_t2star_D1)/sqrt(size(acrossRun_left_ATL_t2star_D1,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_t2star_D1,1),confint,'LineStyle','none','Color','black');
% across run - left ATL - tedana
subplot(6,2,7)
bar(subcode,mean(acrossRun_left_ATL_tedana_D1,1),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ATL - tedana')
hold on
confint = 1.96*(std(acrossRun_left_ATL_tedana_D1)/sqrt(size(acrossRun_left_ATL_tedana_D1,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_tedana_D1,1),confint,'LineStyle','none','Color','black');
% across run - right ATL - t2star
subplot(6,2,9)
bar(subcode,mean(acrossRun_right_ATL_t2star_D1,1),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ATL - t2star')
hold on
confint = 1.96*(std(acrossRun_right_ATL_t2star_D1)/sqrt(size(acrossRun_right_ATL_t2star_D1,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_t2star_D1,1),confint,'LineStyle','none','Color','black');
% across run - right ATL - tedana
subplot(6,2,11)
bar(subcode,mean(acrossRun_right_ATL_tedana_D1,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ATL - tedana')
hold on
confint = 1.96*(std(acrossRun_right_ATL_tedana_D1)/sqrt(size(acrossRun_right_ATL_tedana_D1,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_tedana_D1,1),confint,'LineStyle','none','Color','black');
% averaged - whole brain - t2star
subplot(6,2,2)
bar(subcode,mean(averaged_whole_brain_t2star_D1,1),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Averaged - whole brain - t2star')
hold on
confint = 1.96*(std(averaged_whole_brain_t2star_D1)/sqrt(size(averaged_whole_brain_t2star_D1,1)));
errorbar(1:size(subcode,2),mean(averaged_whole_brain_t2star_D1,1),confint,'LineStyle','none','Color','black');
% averaged - whole brain - tedana
subplot(6,2,4)
bar(subcode,mean(averaged_whole_brain_tedana_D1,1),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - whole brain - tedana')
hold on
confint = 1.96*(std(averaged_whole_brain_tedana_D1)/sqrt(size(averaged_whole_brain_tedana_D1,1)));
errorbar(1:size(subcode,2),mean(averaged_whole_brain_tedana_D1,1),confint,'LineStyle','none','Color','black');
% averaged - left ATL - t2star
subplot(6,2,6)
bar(subcode,mean(averaged_left_ATL_t2star_D1,1),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Averaged - left ATL - t2star')
hold on
confint = 1.96*(std(averaged_left_ATL_t2star_D1)/sqrt(size(averaged_left_ATL_t2star_D1,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_t2star_D1,1),confint,'LineStyle','none','Color','black');
% averaged - left ATL - tedana
subplot(6,2,8)
bar(subcode,mean(averaged_left_ATL_tedana_D1,1),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ATL - tedana')
hold on
confint = 1.96*(std(averaged_left_ATL_tedana_D1)/sqrt(size(averaged_left_ATL_tedana_D1,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_tedana_D1,1),confint,'LineStyle','none','Color','black');
% averaged - right ATL - t2star
subplot(6,2,10)
bar(subcode,mean(averaged_right_ATL_t2star_D1,1),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ATL - t2star')
hold on
confint = 1.96*(std(averaged_right_ATL_t2star_D1)/sqrt(size(averaged_right_ATL_t2star_D1,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_t2star_D1,1),confint,'LineStyle','none','Color','black');
% averaged - right ATL - tedana
subplot(6,2,12)
bar(subcode,mean(averaged_right_ATL_tedana_D1,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ATL - tedana')
hold on
confint = 1.96*(std(averaged_right_ATL_tedana_D1)/sqrt(size(averaged_right_ATL_tedana_D1,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_tedana_D1,1),confint,'LineStyle','none','Color','black');
sgtitle('Linear results - dimension 1')

% second figure shows dimension 2 for every participants
figure;
% across run - whole brain - t2star
subplot(6,2,1)
bar(subcode,mean(acrossRun_whole_brain_t2star_D2,1),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Across run - whole brain - t2star')
hold on
confint = 1.96*(std(acrossRun_whole_brain_t2star_D2)/sqrt(size(acrossRun_whole_brain_t2star_D2,1)));
errorbar(1:size(subcode,2),mean(acrossRun_whole_brain_t2star_D2,1),confint,'LineStyle','none','Color','black');
% across run - whole brain - tedana
subplot(6,2,3)
bar(subcode,mean(acrossRun_whole_brain_tedana_D2,1),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - whole brain - tedana')
hold on
confint = 1.96*(std(acrossRun_whole_brain_tedana_D2)/sqrt(size(acrossRun_whole_brain_tedana_D2,1)));
errorbar(1:size(subcode,2),mean(acrossRun_whole_brain_tedana_D2,1),confint,'LineStyle','none','Color','black');
% across run - left ATL - t2star
subplot(6,2,5)
bar(subcode,mean(acrossRun_left_ATL_t2star_D2,1),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Across run - left ATL - t2star')
hold on
confint = 1.96*(std(acrossRun_left_ATL_t2star_D2)/sqrt(size(acrossRun_left_ATL_t2star_D2,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_t2star_D2,1),confint,'LineStyle','none','Color','black');
% across run - left ATL - tedana
subplot(6,2,7)
bar(subcode,mean(acrossRun_left_ATL_tedana_D2,1),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ATL - tedana')
hold on
confint = 1.96*(std(acrossRun_left_ATL_tedana_D2)/sqrt(size(acrossRun_left_ATL_tedana_D2,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_tedana_D2,1),confint,'LineStyle','none','Color','black');
% across run - right ATL - t2star
subplot(6,2,9)
bar(subcode,mean(acrossRun_right_ATL_t2star_D2,1),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ATL - t2star')
hold on
confint = 1.96*(std(acrossRun_right_ATL_t2star_D2)/sqrt(size(acrossRun_right_ATL_t2star_D2,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_t2star_D2,1),confint,'LineStyle','none','Color','black');
% across run - right ATL - tedana
subplot(6,2,11)
bar(subcode,mean(acrossRun_right_ATL_tedana_D2,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ATL - tedana')
hold on
confint = 1.96*(std(acrossRun_right_ATL_tedana_D2)/sqrt(size(acrossRun_right_ATL_tedana_D2,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_tedana_D2,1),confint,'LineStyle','none','Color','black');
% averaged - whole brain - t2star
subplot(6,2,2)
bar(subcode,mean(averaged_whole_brain_t2star_D2,1),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Averaged - whole brain - t2star')
hold on
confint = 1.96*(std(averaged_whole_brain_t2star_D2)/sqrt(size(averaged_whole_brain_t2star_D2,1)));
errorbar(1:size(subcode,2),mean(averaged_whole_brain_t2star_D2,1),confint,'LineStyle','none','Color','black');
% averaged - whole brain - tedana
subplot(6,2,4)
bar(subcode,mean(averaged_whole_brain_tedana_D2,1),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - whole brain - tedana')
hold on
confint = 1.96*(std(averaged_whole_brain_tedana_D2)/sqrt(size(averaged_whole_brain_tedana_D2,1)));
errorbar(1:size(subcode,2),mean(averaged_whole_brain_tedana_D2,1),confint,'LineStyle','none','Color','black');
% averaged - left ATL - t2star
subplot(6,2,6)
bar(subcode,mean(averaged_left_ATL_t2star_D2,1),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Averaged - left ATL - t2star')
hold on
confint = 1.96*(std(averaged_left_ATL_t2star_D2)/sqrt(size(averaged_left_ATL_t2star_D2,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_t2star_D2,1),confint,'LineStyle','none','Color','black');
% averaged - left ATL - tedana
subplot(6,2,8)
bar(subcode,mean(averaged_left_ATL_tedana_D2,1),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ATL - tedana')
hold on
confint = 1.96*(std(averaged_left_ATL_tedana_D2)/sqrt(size(averaged_left_ATL_tedana_D2,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_tedana_D2,1),confint,'LineStyle','none','Color','black');
% averaged - right ATL - t2star
subplot(6,2,10)
bar(subcode,mean(averaged_right_ATL_t2star_D2,1),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ATL - t2star')
hold on
confint = 1.96*(std(averaged_right_ATL_t2star_D2)/sqrt(size(averaged_right_ATL_t2star_D2,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_t2star_D2,1),confint,'LineStyle','none','Color','black');
% averaged - right ATL - tedana
subplot(6,2,12)
bar(subcode,mean(averaged_right_ATL_tedana_D2,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ATL - tedana')
hold on
confint = 1.96*(std(averaged_right_ATL_tedana_D2)/sqrt(size(averaged_right_ATL_tedana_D2,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_tedana_D2,1),confint,'LineStyle','none','Color','black');
sgtitle('Linear results - dimension 2')

% third figure shows dimension 3 for every participant
figure;
% across run - whole brain - t2star
subplot(6,2,1)
bar(subcode,mean(acrossRun_whole_brain_t2star_D3,1),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Across run - whole brain - t2star')
hold on
confint = 1.96*(std(acrossRun_whole_brain_t2star_D3)/sqrt(size(acrossRun_whole_brain_t2star_D3,1)));
errorbar(1:size(subcode,2),mean(acrossRun_whole_brain_t2star_D3,1),confint,'LineStyle','none','Color','black');
% across run - whole brain - tedana
subplot(6,2,3)
bar(subcode,mean(acrossRun_whole_brain_tedana_D3,1),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - whole brain - tedana')
hold on
confint = 1.96*(std(acrossRun_whole_brain_tedana_D3)/sqrt(size(acrossRun_whole_brain_tedana_D3,1)));
errorbar(1:size(subcode,2),mean(acrossRun_whole_brain_tedana_D3,1),confint,'LineStyle','none','Color','black');
% across run - left ATL - t2star
subplot(6,2,5)
bar(subcode,mean(acrossRun_left_ATL_t2star_D3,1),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Across run - left ATL - t2star')
hold on
confint = 1.96*(std(acrossRun_left_ATL_t2star_D3)/sqrt(size(acrossRun_left_ATL_t2star_D3,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_t2star_D3,1),confint,'LineStyle','none','Color','black');
% across run - left ATL - tedana
subplot(6,2,7)
bar(subcode,mean(acrossRun_left_ATL_tedana_D3,1),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ATL - tedana')
hold on
confint = 1.96*(std(acrossRun_left_ATL_tedana_D3)/sqrt(size(acrossRun_left_ATL_tedana_D3,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_tedana_D3,1),confint,'LineStyle','none','Color','black');
% across run - right ATL - t2star
subplot(6,2,9)
bar(subcode,mean(acrossRun_right_ATL_t2star_D3,1),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ATL - t2star')
hold on
confint = 1.96*(std(acrossRun_right_ATL_t2star_D3)/sqrt(size(acrossRun_right_ATL_t2star_D3,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_t2star_D3,1),confint,'LineStyle','none','Color','black');
% across run - right ATL - tedana
subplot(6,2,11)
bar(subcode,mean(acrossRun_right_ATL_tedana_D3,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ATL - tedana')
hold on
confint = 1.96*(std(acrossRun_right_ATL_tedana_D3)/sqrt(size(acrossRun_right_ATL_tedana_D3,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_tedana_D3,1),confint,'LineStyle','none','Color','black');
% averaged - whole brain - t2star
subplot(6,2,2)
bar(subcode,mean(averaged_whole_brain_t2star_D3,1),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Averaged - whole brain - t2star')
hold on
confint = 1.96*(std(averaged_whole_brain_t2star_D3)/sqrt(size(averaged_whole_brain_t2star_D3,1)));
errorbar(1:size(subcode,2),mean(averaged_whole_brain_t2star_D3,1),confint,'LineStyle','none','Color','black');
% averaged - whole brain - tedana
subplot(6,2,4)
bar(subcode,mean(averaged_whole_brain_tedana_D3,1),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - whole brain - tedana')
hold on
confint = 1.96*(std(averaged_whole_brain_tedana_D3)/sqrt(size(averaged_whole_brain_tedana_D3,1)));
errorbar(1:size(subcode,2),mean(averaged_whole_brain_tedana_D3,1),confint,'LineStyle','none','Color','black');
% averaged - left ATL - t2star
subplot(6,2,6)
bar(subcode,mean(averaged_left_ATL_t2star_D3,1),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Averaged - left ATL - t2star')
hold on
confint = 1.96*(std(averaged_left_ATL_t2star_D3)/sqrt(size(averaged_left_ATL_t2star_D3,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_t2star_D3,1),confint,'LineStyle','none','Color','black');
% averaged - left ATL - tedana
subplot(6,2,8)
bar(subcode,mean(averaged_left_ATL_tedana_D3,1),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ATL - tedana')
hold on
confint = 1.96*(std(averaged_left_ATL_tedana_D3)/sqrt(size(averaged_left_ATL_tedana_D3,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_tedana_D3,1),confint,'LineStyle','none','Color','black');
% averaged - right ATL - t2star
subplot(6,2,10)
bar(subcode,mean(averaged_right_ATL_t2star_D3,1),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ATL - t2star')
hold on
confint = 1.96*(std(averaged_right_ATL_t2star_D3)/sqrt(size(averaged_right_ATL_t2star_D3,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_t2star_D3,1),confint,'LineStyle','none','Color','black');
% averaged - right ATL - tedana
subplot(6,2,12)
bar(subcode,mean(averaged_right_ATL_tedana_D3,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ATL - tedana')
hold on
confint = 1.96*(std(averaged_right_ATL_tedana_D3)/sqrt(size(averaged_right_ATL_tedana_D3,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_tedana_D3,1),confint,'LineStyle','none','Color','black');
sgtitle('Linear results - dimension 3')

% fourth figure shows group averages
group = zeros(size(subcode,2),36);
group(:,1) = mean(acrossRun_whole_brain_t2star_D1)';
group(:,2) = mean(acrossRun_whole_brain_t2star_D2)';
group(:,3) = mean(acrossRun_whole_brain_t2star_D3)';
group(:,4) = mean(acrossRun_whole_brain_tedana_D1)';
group(:,5) = mean(acrossRun_whole_brain_tedana_D2)';
group(:,6) = mean(acrossRun_whole_brain_tedana_D3)';
group(:,7) = mean(acrossRun_left_ATL_t2star_D1)';
group(:,8) = mean(acrossRun_left_ATL_t2star_D2)';
group(:,9) = mean(acrossRun_left_ATL_t2star_D3)';
group(:,10) = mean(acrossRun_left_ATL_tedana_D1)';
group(:,11) = mean(acrossRun_left_ATL_tedana_D2)';
group(:,12) = mean(acrossRun_left_ATL_tedana_D3)';
group(:,13) = mean(acrossRun_right_ATL_t2star_D1)';
group(:,14) = mean(acrossRun_right_ATL_t2star_D2)';
group(:,15) = mean(acrossRun_right_ATL_t2star_D3)';
group(:,16) = mean(acrossRun_right_ATL_tedana_D1)';
group(:,17) = mean(acrossRun_right_ATL_tedana_D2)';
group(:,18) = mean(acrossRun_right_ATL_tedana_D3)';
group(:,19) = mean(averaged_whole_brain_t2star_D1)';
group(:,20) = mean(averaged_whole_brain_t2star_D2)';
group(:,21) = mean(averaged_whole_brain_t2star_D3)';
group(:,22) = mean(averaged_whole_brain_tedana_D1)';
group(:,23) = mean(averaged_whole_brain_tedana_D2)';
group(:,24) = mean(averaged_whole_brain_tedana_D3)';
group(:,25) = mean(averaged_left_ATL_t2star_D1)';
group(:,26) = mean(averaged_left_ATL_t2star_D2)';
group(:,27) = mean(averaged_left_ATL_t2star_D3)';
group(:,28) = mean(averaged_left_ATL_tedana_D1)';
group(:,29) = mean(averaged_left_ATL_tedana_D2)';
group(:,30) = mean(averaged_left_ATL_tedana_D3)';
group(:,31) = mean(averaged_right_ATL_t2star_D1)';
group(:,32) = mean(averaged_right_ATL_t2star_D2)';
group(:,33) = mean(averaged_right_ATL_t2star_D3)';
group(:,34) = mean(averaged_right_ATL_tedana_D1)';
group(:,35) = mean(averaged_right_ATL_tedana_D2)';
group(:,36) = mean(averaged_right_ATL_tedana_D3)';
colnames = {'AR-WB-t2star-D1','AR-WB-t2star-D2','AR-WB-t2star-D3','AR-WB-tedana-D1','AR-WB-tedana-D2','AR-WB-tedana-D3','AR-lATL-t2star-D1','AR-lATL-t2star-D2','AR-lATL-t2star-D3','AR-lATL-tedana-D1','AR-lATL-tedana-D2','AR-lATL-tedana-D3','AR-rATL-t2star-D1','AR-rATL-t2star-D2','AR-rATL-t2star-D3','AR-rATL-tedana-D1','AR-rATL-tedana-D2','AR-rATL-tedana-D3','av-WB-t2star-D1','av-WB-t2star-D2','av-WB-t2star-D3','av-WB-tedana-D1','av-WB-tedana-D2','av-WB-tedana-D3','av-lATL-t2star-D1','av-lATL-t2star-D2','av-lATL-t2star-D3','av-lATL-tedana-D1','av-lATL-tedana-D2','av-lATL-tedana-D3','av-rATL-t2star-D1','av-rATL-t2star-D2','av-rATL-t2star-D3','av-rATL-tedana-D1','av-rATL-tedana-D2','av-rATL-tedana-D3'};

figure;
b = bar(colnames,mean(group),'FaceColor','flat');
% normalised RGB triplets
colours = [119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189;119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189]/255;
b.CData = colours;
ylim([-0.1,1])
title('Group results')
hold on
confint = 1.96*(std(group)/sqrt(size(group,1)));
errorbar(1:size(group,2),mean(group,1),confint,'LineStyle','none','Color','black');

% fifth figure shows results for anterior and posterior ROIs for every
% participant on dimension 1
figure;
% across run - left ROI (anterior) - t2star
subplot(8,2,1)
bar(subcode,mean(acrossRun_left_ATL_ant_t2star_D1,1),'FaceColor','#D95319')
ylim([-0.1,1])
title('Across run - left ROI (anterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_left_ATL_ant_t2star_D1)/sqrt(size(acrossRun_left_ATL_ant_t2star_D1,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_ant_t2star_D1,1),confint,'LineStyle','none','Color','black');
% across run - left ROI (anterior) - tedana
subplot(8,2,3)
bar(subcode,mean(acrossRun_left_ATL_ant_tedana_D1,1),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_left_ATL_ant_tedana_D1)/sqrt(size(acrossRun_left_ATL_ant_tedana_D1,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_ant_tedana_D1,1),confint,'LineStyle','none','Color','black');
% across run - left ROI (posterior) - t2star
subplot(8,2,5)
bar(subcode,mean(acrossRun_left_ATL_pos_t2star_D1,1),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Across run - left ROI (posterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_left_ATL_pos_t2star_D1)/sqrt(size(acrossRun_left_ATL_pos_t2star_D1,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_pos_t2star_D1,1),confint,'LineStyle','none','Color','black');
% across run - left ROI (posterior) - tedana
subplot(8,2,7)
bar(subcode,mean(acrossRun_left_ATL_pos_tedana_D1,1),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_left_ATL_pos_tedana_D1)/sqrt(size(acrossRun_left_ATL_pos_tedana_D1,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_pos_tedana_D1,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (anterior) - t2star
subplot(8,2,9)
bar(subcode,mean(acrossRun_right_ATL_ant_t2star_D1,1),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ROI (anterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_right_ATL_ant_t2star_D1)/sqrt(size(acrossRun_right_ATL_ant_t2star_D1,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_ant_t2star_D1,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (anterior) - tedana
subplot(8,2,11)
bar(subcode,mean(acrossRun_right_ATL_ant_tedana_D1,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_right_ATL_ant_tedana_D1)/sqrt(size(acrossRun_right_ATL_ant_tedana_D1,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_ant_tedana_D1,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (posterior) - t2star
subplot(8,2,13)
bar(subcode,mean(acrossRun_right_ATL_pos_t2star_D1,1),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Across run - right ROI (posterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_right_ATL_pos_t2star_D1)/sqrt(size(acrossRun_right_ATL_pos_t2star_D1,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_pos_t2star_D1,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (posterior) - tedana
subplot(8,2,15)
bar(subcode,mean(acrossRun_right_ATL_pos_tedana_D1,1),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_right_ATL_pos_tedana_D1)/sqrt(size(acrossRun_right_ATL_pos_tedana_D1,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_pos_tedana_D1,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (anterior) - t2star
subplot(8,2,2)
bar(subcode,mean(averaged_left_ATL_ant_t2star_D1,1),'FaceColor','#D95319')
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - t2star')
hold on
confint = 1.96*(std(averaged_left_ATL_ant_t2star_D1)/sqrt(size(averaged_left_ATL_ant_t2star_D1,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_ant_t2star_D1,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (anterior) - tedana
subplot(8,2,4)
bar(subcode,mean(averaged_left_ATL_ant_tedana_D1,1),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - tedana')
hold on
confint = 1.96*(std(averaged_left_ATL_ant_tedana_D1)/sqrt(size(averaged_left_ATL_ant_tedana_D1,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_ant_tedana_D1,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (posterior) - t2star
subplot(8,2,6)
bar(subcode,mean(averaged_left_ATL_pos_t2star_D1,1),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - t2star')
hold on
confint = 1.96*(std(averaged_left_ATL_pos_t2star_D1)/sqrt(size(averaged_left_ATL_pos_t2star_D1,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_pos_t2star_D1,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (posterior) - tedana
subplot(8,2,8)
bar(subcode,mean(averaged_left_ATL_pos_tedana_D1,1),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - tedana')
hold on
confint = 1.96*(std(averaged_left_ATL_pos_tedana_D1)/sqrt(size(averaged_left_ATL_pos_tedana_D1,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_pos_tedana_D1,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (anterior) - t2star
subplot(8,2,10)
bar(subcode,mean(averaged_right_ATL_ant_t2star_D1,1),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - t2star')
hold on
confint = 1.96*(std(averaged_right_ATL_ant_t2star_D1)/sqrt(size(averaged_right_ATL_ant_t2star_D1,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_ant_t2star_D1,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (anterior) - tedana
subplot(8,2,12)
bar(subcode,mean(averaged_right_ATL_ant_tedana_D1,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(averaged_right_ATL_ant_tedana_D1)/sqrt(size(averaged_right_ATL_ant_tedana_D1,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_ant_tedana_D1,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (posterior) - t2star
subplot(8,2,14)
bar(subcode,mean(averaged_right_ATL_pos_t2star_D1,1),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Averaged - right ROI (posterior) - t2star')
hold on
confint = 1.96*(std(averaged_right_ATL_pos_t2star_D1)/sqrt(size(averaged_right_ATL_pos_t2star_D1,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_pos_t2star_D1,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (posterior) - tedana
subplot(8,2,16)
bar(subcode,mean(averaged_right_ATL_pos_tedana_D1,1),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(averaged_right_ATL_pos_tedana_D1)/sqrt(size(averaged_right_ATL_pos_tedana_D1,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_pos_tedana_D1,1),confint,'LineStyle','none','Color','black');
sgtitle('Dimension 1 - anterior and posterior halves of ROI')

% sixth figure shows results for anterior and posterior ROIs for every
% participant on dimension 2
figure;
% across run - left ROI (anterior) - t2star
subplot(8,2,1)
bar(subcode,mean(acrossRun_left_ATL_ant_t2star_D2,1),'FaceColor','#D95319')
ylim([-0.1,1])
title('Across run - left ROI (anterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_left_ATL_ant_t2star_D2)/sqrt(size(acrossRun_left_ATL_ant_t2star_D2,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_ant_t2star_D2,1),confint,'LineStyle','none','Color','black');
% across run - left ROI (anterior) - tedana
subplot(8,2,3)
bar(subcode,mean(acrossRun_left_ATL_ant_tedana_D2,1),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_left_ATL_ant_tedana_D2)/sqrt(size(acrossRun_left_ATL_ant_tedana_D2,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_ant_tedana_D2,1),confint,'LineStyle','none','Color','black');
% across run - left ROI (posterior) - t2star
subplot(8,2,5)
bar(subcode,mean(acrossRun_left_ATL_pos_t2star_D2,1),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Across run - left ROI (posterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_left_ATL_pos_t2star_D2)/sqrt(size(acrossRun_left_ATL_pos_t2star_D2,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_pos_t2star_D2,1),confint,'LineStyle','none','Color','black');
% across run - left ROI (posterior) - tedana
subplot(8,2,7)
bar(subcode,mean(acrossRun_left_ATL_pos_tedana_D2,1),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_left_ATL_pos_tedana_D2)/sqrt(size(acrossRun_left_ATL_pos_tedana_D2,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_pos_tedana_D2,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (anterior) - t2star
subplot(8,2,9)
bar(subcode,mean(acrossRun_right_ATL_ant_t2star_D2,1),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ROI (anterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_right_ATL_ant_t2star_D2)/sqrt(size(acrossRun_right_ATL_ant_t2star_D2,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_ant_t2star_D2,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (anterior) - tedana
subplot(8,2,11)
bar(subcode,mean(acrossRun_right_ATL_ant_tedana_D2,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_right_ATL_ant_tedana_D2)/sqrt(size(acrossRun_right_ATL_ant_tedana_D2,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_ant_tedana_D2,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (posterior) - t2star
subplot(8,2,13)
bar(subcode,mean(acrossRun_right_ATL_pos_t2star_D2,1),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Across run - right ROI (posterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_right_ATL_pos_t2star_D2)/sqrt(size(acrossRun_right_ATL_pos_t2star_D2,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_pos_t2star_D2,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (posterior) - tedana
subplot(8,2,15)
bar(subcode,mean(acrossRun_right_ATL_pos_tedana_D2,1),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_right_ATL_pos_tedana_D2)/sqrt(size(acrossRun_right_ATL_pos_tedana_D2,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_pos_tedana_D2,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (anterior) - t2star
subplot(8,2,2)
bar(subcode,mean(averaged_left_ATL_ant_t2star_D2,1),'FaceColor','#D95319')
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - t2star')
hold on
confint = 1.96*(std(averaged_left_ATL_ant_t2star_D2)/sqrt(size(averaged_left_ATL_ant_t2star_D2,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_ant_t2star_D2,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (anterior) - tedana
subplot(8,2,4)
bar(subcode,mean(averaged_left_ATL_ant_tedana_D2,1),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - tedana')
hold on
confint = 1.96*(std(averaged_left_ATL_ant_tedana_D2)/sqrt(size(averaged_left_ATL_ant_tedana_D2,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_ant_tedana_D2,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (posterior) - t2star
subplot(8,2,6)
bar(subcode,mean(averaged_left_ATL_pos_t2star_D2,1),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - t2star')
hold on
confint = 1.96*(std(averaged_left_ATL_pos_t2star_D2)/sqrt(size(averaged_left_ATL_pos_t2star_D2,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_pos_t2star_D2,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (posterior) - tedana
subplot(8,2,8)
bar(subcode,mean(averaged_left_ATL_pos_tedana_D2,1),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - tedana')
hold on
confint = 1.96*(std(averaged_left_ATL_pos_tedana_D2)/sqrt(size(averaged_left_ATL_pos_tedana_D2,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_pos_tedana_D2,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (anterior) - t2star
subplot(8,2,10)
bar(subcode,mean(averaged_right_ATL_ant_t2star_D2,1),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - t2star')
hold on
confint = 1.96*(std(averaged_right_ATL_ant_t2star_D2)/sqrt(size(averaged_right_ATL_ant_t2star_D2,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_ant_t2star_D2,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (anterior) - tedana
subplot(8,2,12)
bar(subcode,mean(averaged_right_ATL_ant_tedana_D2,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(averaged_right_ATL_ant_tedana_D2)/sqrt(size(averaged_right_ATL_ant_tedana_D2,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_ant_tedana_D2,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (posterior) - t2star
subplot(8,2,14)
bar(subcode,mean(averaged_right_ATL_pos_t2star_D2,1),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Averaged - right ROI (posterior) - t2star')
hold on
confint = 1.96*(std(averaged_right_ATL_pos_t2star_D2)/sqrt(size(averaged_right_ATL_pos_t2star_D2,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_pos_t2star_D2,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (posterior) - tedana
subplot(8,2,16)
bar(subcode,mean(averaged_right_ATL_pos_tedana_D2,1),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(averaged_right_ATL_pos_tedana_D2)/sqrt(size(averaged_right_ATL_pos_tedana_D2,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_pos_tedana_D2,1),confint,'LineStyle','none','Color','black');
sgtitle('Dimension 2 - anterior and posterior halves of ROI')

% seventh figure shows results for anterior and posterior ROIs for every
% participant on dimension 3
figure;
% across run - left ROI (anterior) - t2star
subplot(8,2,1)
bar(subcode,mean(acrossRun_left_ATL_ant_t2star_D3,1),'FaceColor','#D95319')
ylim([-0.1,1])
title('Across run - left ROI (anterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_left_ATL_ant_t2star_D3)/sqrt(size(acrossRun_left_ATL_ant_t2star_D3,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_ant_t2star_D3,1),confint,'LineStyle','none','Color','black');
% across run - left ROI (anterior) - tedana
subplot(8,2,3)
bar(subcode,mean(acrossRun_left_ATL_ant_tedana_D3,1),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_left_ATL_ant_tedana_D3)/sqrt(size(acrossRun_left_ATL_ant_tedana_D3,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_ant_tedana_D3,1),confint,'LineStyle','none','Color','black');
% across run - left ROI (posterior) - t2star
subplot(8,2,5)
bar(subcode,mean(acrossRun_left_ATL_pos_t2star_D3,1),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Across run - left ROI (posterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_left_ATL_pos_t2star_D3)/sqrt(size(acrossRun_left_ATL_pos_t2star_D3,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_pos_t2star_D3,1),confint,'LineStyle','none','Color','black');
% across run - left ROI (posterior) - tedana
subplot(8,2,7)
bar(subcode,mean(acrossRun_left_ATL_pos_tedana_D3,1),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_left_ATL_pos_tedana_D3)/sqrt(size(acrossRun_left_ATL_pos_tedana_D3,1)));
errorbar(1:size(subcode,2),mean(acrossRun_left_ATL_pos_tedana_D3,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (anterior) - t2star
subplot(8,2,9)
bar(subcode,mean(acrossRun_right_ATL_ant_t2star_D3,1),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ROI (anterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_right_ATL_ant_t2star_D3)/sqrt(size(acrossRun_right_ATL_ant_t2star_D3,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_ant_t2star_D3,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (anterior) - tedana
subplot(8,2,11)
bar(subcode,mean(acrossRun_right_ATL_ant_tedana_D3,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_right_ATL_ant_tedana_D3)/sqrt(size(acrossRun_right_ATL_ant_tedana_D3,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_ant_tedana_D3,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (posterior) - t2star
subplot(8,2,13)
bar(subcode,mean(acrossRun_right_ATL_pos_t2star_D3,1),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Across run - right ROI (posterior) - t2star')
hold on
confint = 1.96*(std(acrossRun_right_ATL_pos_t2star_D3)/sqrt(size(acrossRun_right_ATL_pos_t2star_D3,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_pos_t2star_D3,1),confint,'LineStyle','none','Color','black');
% across run - right ROI (posterior) - tedana
subplot(8,2,15)
bar(subcode,mean(acrossRun_right_ATL_pos_tedana_D3,1),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(acrossRun_right_ATL_pos_tedana_D3)/sqrt(size(acrossRun_right_ATL_pos_tedana_D3,1)));
errorbar(1:size(subcode,2),mean(acrossRun_right_ATL_pos_tedana_D3,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (anterior) - t2star
subplot(8,2,2)
bar(subcode,mean(averaged_left_ATL_ant_t2star_D3,1),'FaceColor','#D95319')
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - t2star')
hold on
confint = 1.96*(std(averaged_left_ATL_ant_t2star_D3)/sqrt(size(averaged_left_ATL_ant_t2star_D3,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_ant_t2star_D3,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (anterior) - tedana
subplot(8,2,4)
bar(subcode,mean(averaged_left_ATL_ant_tedana_D3,1),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - tedana')
hold on
confint = 1.96*(std(averaged_left_ATL_ant_tedana_D3)/sqrt(size(averaged_left_ATL_ant_tedana_D3,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_ant_tedana_D3,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (posterior) - t2star
subplot(8,2,6)
bar(subcode,mean(averaged_left_ATL_pos_t2star_D3,1),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - t2star')
hold on
confint = 1.96*(std(averaged_left_ATL_pos_t2star_D3)/sqrt(size(averaged_left_ATL_pos_t2star_D3,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_pos_t2star_D3,1),confint,'LineStyle','none','Color','black');
% averaged - left ROI (posterior) - tedana
subplot(8,2,8)
bar(subcode,mean(averaged_left_ATL_pos_tedana_D3,1),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - tedana')
hold on
confint = 1.96*(std(averaged_left_ATL_pos_tedana_D3)/sqrt(size(averaged_left_ATL_pos_tedana_D3,1)));
errorbar(1:size(subcode,2),mean(averaged_left_ATL_pos_tedana_D3,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (anterior) - t2star
subplot(8,2,10)
bar(subcode,mean(averaged_right_ATL_ant_t2star_D3,1),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - t2star')
hold on
confint = 1.96*(std(averaged_right_ATL_ant_t2star_D3)/sqrt(size(averaged_right_ATL_ant_t2star_D3,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_ant_t2star_D3,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (anterior) - tedana
subplot(8,2,12)
bar(subcode,mean(averaged_right_ATL_ant_tedana_D3,1),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(averaged_right_ATL_ant_tedana_D3)/sqrt(size(averaged_right_ATL_ant_tedana_D3,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_ant_tedana_D3,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (posterior) - t2star
subplot(8,2,14)
bar(subcode,mean(averaged_right_ATL_pos_t2star_D3,1),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Averaged - right ROI (posterior) - t2star')
hold on
confint = 1.96*(std(averaged_right_ATL_pos_t2star_D3)/sqrt(size(averaged_right_ATL_pos_t2star_D3,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_pos_t2star_D3,1),confint,'LineStyle','none','Color','black');
% averaged - right ROI (posterior) - tedana
subplot(8,2,16)
bar(subcode,mean(averaged_right_ATL_pos_tedana_D3,1),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
hold on
confint = 1.96*(std(averaged_right_ATL_pos_tedana_D3)/sqrt(size(averaged_right_ATL_pos_tedana_D3,1)));
errorbar(1:size(subcode,2),mean(averaged_right_ATL_pos_tedana_D3,1),confint,'LineStyle','none','Color','black');
sgtitle('Dimension 3 - anterior and posterior halves of ROI')

% eighth figure shows group results for the above comparisons
groupsubdivide = zeros(size(subcode,2),48);
groupsubdivide(:,1) = mean(acrossRun_left_ATL_ant_t2star_D1)';
groupsubdivide(:,2) = mean(acrossRun_left_ATL_ant_t2star_D2)';
groupsubdivide(:,3) = mean(acrossRun_left_ATL_ant_t2star_D3)';
groupsubdivide(:,4) = mean(acrossRun_left_ATL_ant_tedana_D1)';
groupsubdivide(:,5) = mean(acrossRun_left_ATL_ant_tedana_D2)';
groupsubdivide(:,6) = mean(acrossRun_left_ATL_ant_tedana_D3)';
groupsubdivide(:,7) = mean(acrossRun_left_ATL_pos_t2star_D1)';
groupsubdivide(:,8) = mean(acrossRun_left_ATL_pos_t2star_D2)';
groupsubdivide(:,9) = mean(acrossRun_left_ATL_pos_t2star_D3)';
groupsubdivide(:,10) = mean(acrossRun_left_ATL_pos_tedana_D1)';
groupsubdivide(:,11) = mean(acrossRun_left_ATL_pos_tedana_D2)';
groupsubdivide(:,12) = mean(acrossRun_left_ATL_pos_tedana_D3)';
groupsubdivide(:,13) = mean(acrossRun_right_ATL_ant_t2star_D1)';
groupsubdivide(:,14) = mean(acrossRun_right_ATL_ant_t2star_D2)';
groupsubdivide(:,15) = mean(acrossRun_right_ATL_ant_t2star_D3)';
groupsubdivide(:,16) = mean(acrossRun_right_ATL_ant_tedana_D1)';
groupsubdivide(:,17) = mean(acrossRun_right_ATL_ant_tedana_D2)';
groupsubdivide(:,18) = mean(acrossRun_right_ATL_ant_tedana_D3)';
groupsubdivide(:,19) = mean(acrossRun_right_ATL_pos_t2star_D1)';
groupsubdivide(:,20) = mean(acrossRun_right_ATL_pos_t2star_D2)';
groupsubdivide(:,21) = mean(acrossRun_right_ATL_pos_t2star_D3)';
groupsubdivide(:,22) = mean(acrossRun_right_ATL_pos_tedana_D1)';
groupsubdivide(:,23) = mean(acrossRun_right_ATL_pos_tedana_D2)';
groupsubdivide(:,24) = mean(acrossRun_right_ATL_pos_tedana_D3)';
groupsubdivide(:,25) = mean(averaged_left_ATL_ant_t2star_D1)';
groupsubdivide(:,26) = mean(averaged_left_ATL_ant_t2star_D2)';
groupsubdivide(:,27) = mean(averaged_left_ATL_ant_t2star_D3)';
groupsubdivide(:,28) = mean(averaged_left_ATL_ant_tedana_D1)';
groupsubdivide(:,29) = mean(averaged_left_ATL_ant_tedana_D2)';
groupsubdivide(:,30) = mean(averaged_left_ATL_ant_tedana_D3)';
groupsubdivide(:,31) = mean(averaged_left_ATL_pos_t2star_D1)';
groupsubdivide(:,32) = mean(averaged_left_ATL_pos_t2star_D2)';
groupsubdivide(:,33) = mean(averaged_left_ATL_pos_t2star_D3)';
groupsubdivide(:,34) = mean(averaged_left_ATL_pos_tedana_D1)';
groupsubdivide(:,35) = mean(averaged_left_ATL_pos_tedana_D2)';
groupsubdivide(:,36) = mean(averaged_left_ATL_pos_tedana_D3)';
groupsubdivide(:,37) = mean(averaged_right_ATL_ant_t2star_D1)';
groupsubdivide(:,38) = mean(averaged_right_ATL_ant_t2star_D2)';
groupsubdivide(:,39) = mean(averaged_right_ATL_ant_t2star_D3)';
groupsubdivide(:,40) = mean(averaged_right_ATL_ant_tedana_D1)';
groupsubdivide(:,41) = mean(averaged_right_ATL_ant_tedana_D2)';
groupsubdivide(:,42) = mean(averaged_right_ATL_ant_tedana_D3)';
groupsubdivide(:,43) = mean(averaged_right_ATL_pos_t2star_D1)';
groupsubdivide(:,44) = mean(averaged_right_ATL_pos_t2star_D2)';
groupsubdivide(:,45) = mean(averaged_right_ATL_pos_t2star_D3)';
groupsubdivide(:,46) = mean(averaged_right_ATL_pos_tedana_D1)';
groupsubdivide(:,47) = mean(averaged_right_ATL_pos_tedana_D2)';
groupsubdivide(:,48) = mean(averaged_right_ATL_pos_tedana_D3)';
colnames = {'AR-lATL-ant-t2star-D1','AR-lATL-ant-t2star-D2','AR-lATL-ant-t2star-D3','AR-lATL-ant-tedana-D1','AR-lATL-ant-tedana-D2','AR-lATL-ant-tedana-D3','AR-lATL-pos-t2star-D1','AR-lATL-pos-t2star-D2','AR-lATL-pos-t2star-D3','AR-lATL-pos-tedana-D1','AR-lATL-pos-tedana-D2','AR-lATL-pos-tedana-D3','AR-rATL-ant-t2star-D1','AR-rATL-ant-t2star-D2','AR-rATL-ant-t2star-D3','AR-rATL-ant-tedana-D1','AR-rATL-ant-tedana-D2','AR-rATL-ant-tedana-D3','AR-rATL-pos-t2star-D1','AR-rATL-pos-t2star-D2','AR-rATL-pos-t2star-D3','AR-rATL-pos-tedana-D1','AR-rATL-pos-tedana-D2','AR-rATL-pos-tedana-D3','av-lATL-ant-t2star-D1','av-lATL-ant-t2star-D2','av-lATL-ant-t2star-D3','av-lATL-ant-tedana-D1','av-lATL-ant-tedana-D2','av-lATL-ant-tedana-D3','av-lATL-pos-t2star-D1','av-lATL-pos-t2star-D2','av-lATL-pos-t2star-D3','av-lATL-pos-tedana-D1','av-lATL-pos-tedana-D2','av-lATL-pos-tedana-D3','av-rATL-ant-t2star-D1','av-rATL-ant-t2star-D2','av-rATL-ant-t2star-D3','av-rATL-ant-tedana-D1','av-rATL-ant-tedana-D2','av-rATL-ant-tedana-D3','av-rATL-pos-t2star-D1','av-rATL-pos-t2star-D2','av-rATL-pos-t2star-D3','av-rATL-pos-tedana-D1','av-rATL-pos-tedana-D2','av-rATL-pos-tedana-D3',};

figure;
b = bar(colnames,mean(groupsubdivide),'FaceColor','flat');
% normalised RGB triplets
colours = [0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330];
b.CData = colours;
ylim([-0.1,1])
title('Group results')
hold on
confint = 1.96*(std(groupsubdivide)/sqrt(size(groupsubdivide,1)));
errorbar(1:size(groupsubdivide,2),mean(groupsubdivide,1),confint,'LineStyle','none','Color','black');

%% correlations with predicted coordinates

% decompose similarity matrix into 3 singular values
load('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/dilkina_norms.mat');
U = embed_similarity_matrix(dilkina_norms,3);

% across-run decoding predicts the same coordinates 4 times. Average these
% predictions
coords_acrossRun_whole_brain_t2star_D1 = mean(cat(3,coords_acrossRun_whole_brain_t2star_D1(1:100,:),coords_acrossRun_whole_brain_t2star_D1(101:200,:),coords_acrossRun_whole_brain_t2star_D1(201:300,:),coords_acrossRun_whole_brain_t2star_D1(301:400,:)),3);
coords_acrossRun_whole_brain_t2star_D2 = mean(cat(3,coords_acrossRun_whole_brain_t2star_D2(1:100,:),coords_acrossRun_whole_brain_t2star_D2(101:200,:),coords_acrossRun_whole_brain_t2star_D2(201:300,:),coords_acrossRun_whole_brain_t2star_D2(301:400,:)),3);
coords_acrossRun_whole_brain_t2star_D3 = mean(cat(3,coords_acrossRun_whole_brain_t2star_D3(1:100,:),coords_acrossRun_whole_brain_t2star_D3(101:200,:),coords_acrossRun_whole_brain_t2star_D3(201:300,:),coords_acrossRun_whole_brain_t2star_D3(301:400,:)),3);
coords_acrossRun_left_ATL_t2star_D1 = mean(cat(3,coords_acrossRun_left_ATL_t2star_D1(1:100,:),coords_acrossRun_left_ATL_t2star_D1(101:200,:),coords_acrossRun_left_ATL_t2star_D1(201:300,:),coords_acrossRun_left_ATL_t2star_D1(301:400,:)),3);
coords_acrossRun_left_ATL_t2star_D2 = mean(cat(3,coords_acrossRun_left_ATL_t2star_D2(1:100,:),coords_acrossRun_left_ATL_t2star_D2(101:200,:),coords_acrossRun_left_ATL_t2star_D2(201:300,:),coords_acrossRun_left_ATL_t2star_D2(301:400,:)),3);
coords_acrossRun_left_ATL_t2star_D3 = mean(cat(3,coords_acrossRun_left_ATL_t2star_D3(1:100,:),coords_acrossRun_left_ATL_t2star_D3(101:200,:),coords_acrossRun_left_ATL_t2star_D3(201:300,:),coords_acrossRun_left_ATL_t2star_D3(301:400,:)),3);
coords_acrossRun_right_ATL_t2star_D1 = mean(cat(3,coords_acrossRun_right_ATL_t2star_D1(1:100,:),coords_acrossRun_right_ATL_t2star_D1(101:200,:),coords_acrossRun_right_ATL_t2star_D1(201:300,:),coords_acrossRun_right_ATL_t2star_D1(301:400,:)),3);
coords_acrossRun_right_ATL_t2star_D2 = mean(cat(3,coords_acrossRun_right_ATL_t2star_D2(1:100,:),coords_acrossRun_right_ATL_t2star_D2(101:200,:),coords_acrossRun_right_ATL_t2star_D2(201:300,:),coords_acrossRun_right_ATL_t2star_D2(301:400,:)),3);
coords_acrossRun_right_ATL_t2star_D3 = mean(cat(3,coords_acrossRun_right_ATL_t2star_D3(1:100,:),coords_acrossRun_right_ATL_t2star_D3(101:200,:),coords_acrossRun_right_ATL_t2star_D3(201:300,:),coords_acrossRun_right_ATL_t2star_D3(301:400,:)),3);
coords_acrossRun_left_ATL_ant_t2star_D1 = mean(cat(3,coords_acrossRun_left_ATL_ant_t2star_D1(1:100,:),coords_acrossRun_left_ATL_ant_t2star_D1(101:200,:),coords_acrossRun_left_ATL_ant_t2star_D1(201:300,:),coords_acrossRun_left_ATL_ant_t2star_D1(301:400,:)),3);
coords_acrossRun_left_ATL_ant_t2star_D2 = mean(cat(3,coords_acrossRun_left_ATL_ant_t2star_D2(1:100,:),coords_acrossRun_left_ATL_ant_t2star_D2(101:200,:),coords_acrossRun_left_ATL_ant_t2star_D2(201:300,:),coords_acrossRun_left_ATL_ant_t2star_D2(301:400,:)),3);
coords_acrossRun_left_ATL_ant_t2star_D3 = mean(cat(3,coords_acrossRun_left_ATL_ant_t2star_D3(1:100,:),coords_acrossRun_left_ATL_ant_t2star_D3(101:200,:),coords_acrossRun_left_ATL_ant_t2star_D3(201:300,:),coords_acrossRun_left_ATL_ant_t2star_D3(301:400,:)),3);
coords_acrossRun_left_ATL_pos_t2star_D1 = mean(cat(3,coords_acrossRun_left_ATL_pos_t2star_D1(1:100,:),coords_acrossRun_left_ATL_pos_t2star_D1(101:200,:),coords_acrossRun_left_ATL_pos_t2star_D1(201:300,:),coords_acrossRun_left_ATL_pos_t2star_D1(301:400,:)),3);
coords_acrossRun_left_ATL_pos_t2star_D2 = mean(cat(3,coords_acrossRun_left_ATL_pos_t2star_D2(1:100,:),coords_acrossRun_left_ATL_pos_t2star_D2(101:200,:),coords_acrossRun_left_ATL_pos_t2star_D2(201:300,:),coords_acrossRun_left_ATL_pos_t2star_D2(301:400,:)),3);
coords_acrossRun_left_ATL_pos_t2star_D3 = mean(cat(3,coords_acrossRun_left_ATL_pos_t2star_D3(1:100,:),coords_acrossRun_left_ATL_pos_t2star_D3(101:200,:),coords_acrossRun_left_ATL_pos_t2star_D3(201:300,:),coords_acrossRun_left_ATL_pos_t2star_D3(301:400,:)),3);
coords_acrossRun_right_ATL_ant_t2star_D1 = mean(cat(3,coords_acrossRun_right_ATL_ant_t2star_D1(1:100,:),coords_acrossRun_right_ATL_ant_t2star_D1(101:200,:),coords_acrossRun_right_ATL_ant_t2star_D1(201:300,:),coords_acrossRun_right_ATL_ant_t2star_D1(301:400,:)),3);
coords_acrossRun_right_ATL_ant_t2star_D2 = mean(cat(3,coords_acrossRun_right_ATL_ant_t2star_D2(1:100,:),coords_acrossRun_right_ATL_ant_t2star_D2(101:200,:),coords_acrossRun_right_ATL_ant_t2star_D2(201:300,:),coords_acrossRun_right_ATL_ant_t2star_D2(301:400,:)),3);
coords_acrossRun_right_ATL_ant_t2star_D3 = mean(cat(3,coords_acrossRun_right_ATL_ant_t2star_D3(1:100,:),coords_acrossRun_right_ATL_ant_t2star_D3(101:200,:),coords_acrossRun_right_ATL_ant_t2star_D3(201:300,:),coords_acrossRun_right_ATL_ant_t2star_D3(301:400,:)),3);
coords_acrossRun_right_ATL_pos_t2star_D1 = mean(cat(3,coords_acrossRun_right_ATL_pos_t2star_D1(1:100,:),coords_acrossRun_right_ATL_pos_t2star_D1(101:200,:),coords_acrossRun_right_ATL_pos_t2star_D1(201:300,:),coords_acrossRun_right_ATL_pos_t2star_D1(301:400,:)),3);
coords_acrossRun_right_ATL_pos_t2star_D2 = mean(cat(3,coords_acrossRun_right_ATL_pos_t2star_D2(1:100,:),coords_acrossRun_right_ATL_pos_t2star_D2(101:200,:),coords_acrossRun_right_ATL_pos_t2star_D2(201:300,:),coords_acrossRun_right_ATL_pos_t2star_D2(301:400,:)),3);
coords_acrossRun_right_ATL_pos_t2star_D3 = mean(cat(3,coords_acrossRun_right_ATL_pos_t2star_D3(1:100,:),coords_acrossRun_right_ATL_pos_t2star_D3(101:200,:),coords_acrossRun_right_ATL_pos_t2star_D3(201:300,:),coords_acrossRun_right_ATL_pos_t2star_D3(301:400,:)),3);
coords_acrossRun_whole_brain_tedana_D1 = mean(cat(3,coords_acrossRun_whole_brain_tedana_D1(1:100,:),coords_acrossRun_whole_brain_tedana_D1(101:200,:),coords_acrossRun_whole_brain_tedana_D1(201:300,:),coords_acrossRun_whole_brain_tedana_D1(301:400,:)),3);
coords_acrossRun_whole_brain_tedana_D2 = mean(cat(3,coords_acrossRun_whole_brain_tedana_D2(1:100,:),coords_acrossRun_whole_brain_tedana_D2(101:200,:),coords_acrossRun_whole_brain_tedana_D2(201:300,:),coords_acrossRun_whole_brain_tedana_D2(301:400,:)),3);
coords_acrossRun_whole_brain_tedana_D3 = mean(cat(3,coords_acrossRun_whole_brain_tedana_D3(1:100,:),coords_acrossRun_whole_brain_tedana_D3(101:200,:),coords_acrossRun_whole_brain_tedana_D3(201:300,:),coords_acrossRun_whole_brain_tedana_D3(301:400,:)),3);
coords_acrossRun_left_ATL_tedana_D1 = mean(cat(3,coords_acrossRun_left_ATL_tedana_D1(1:100,:),coords_acrossRun_left_ATL_tedana_D1(101:200,:),coords_acrossRun_left_ATL_tedana_D1(201:300,:),coords_acrossRun_left_ATL_tedana_D1(301:400,:)),3);
coords_acrossRun_left_ATL_tedana_D2 = mean(cat(3,coords_acrossRun_left_ATL_tedana_D2(1:100,:),coords_acrossRun_left_ATL_tedana_D2(101:200,:),coords_acrossRun_left_ATL_tedana_D2(201:300,:),coords_acrossRun_left_ATL_tedana_D2(301:400,:)),3);
coords_acrossRun_left_ATL_tedana_D3 = mean(cat(3,coords_acrossRun_left_ATL_tedana_D3(1:100,:),coords_acrossRun_left_ATL_tedana_D3(101:200,:),coords_acrossRun_left_ATL_tedana_D3(201:300,:),coords_acrossRun_left_ATL_tedana_D3(301:400,:)),3);
coords_acrossRun_right_ATL_tedana_D1 = mean(cat(3,coords_acrossRun_right_ATL_tedana_D1(1:100,:),coords_acrossRun_right_ATL_tedana_D1(101:200,:),coords_acrossRun_right_ATL_tedana_D1(201:300,:),coords_acrossRun_right_ATL_tedana_D1(301:400,:)),3);
coords_acrossRun_right_ATL_tedana_D2 = mean(cat(3,coords_acrossRun_right_ATL_tedana_D2(1:100,:),coords_acrossRun_right_ATL_tedana_D2(101:200,:),coords_acrossRun_right_ATL_tedana_D2(201:300,:),coords_acrossRun_right_ATL_tedana_D2(301:400,:)),3);
coords_acrossRun_right_ATL_tedana_D3 = mean(cat(3,coords_acrossRun_right_ATL_tedana_D3(1:100,:),coords_acrossRun_right_ATL_tedana_D3(101:200,:),coords_acrossRun_right_ATL_tedana_D3(201:300,:),coords_acrossRun_right_ATL_tedana_D3(301:400,:)),3);
coords_acrossRun_left_ATL_ant_tedana_D1 = mean(cat(3,coords_acrossRun_left_ATL_ant_tedana_D1(1:100,:),coords_acrossRun_left_ATL_ant_tedana_D1(101:200,:),coords_acrossRun_left_ATL_ant_tedana_D1(201:300,:),coords_acrossRun_left_ATL_ant_tedana_D1(301:400,:)),3);
coords_acrossRun_left_ATL_ant_tedana_D2 = mean(cat(3,coords_acrossRun_left_ATL_ant_tedana_D2(1:100,:),coords_acrossRun_left_ATL_ant_tedana_D2(101:200,:),coords_acrossRun_left_ATL_ant_tedana_D2(201:300,:),coords_acrossRun_left_ATL_ant_tedana_D2(301:400,:)),3);
coords_acrossRun_left_ATL_ant_tedana_D3 = mean(cat(3,coords_acrossRun_left_ATL_ant_tedana_D3(1:100,:),coords_acrossRun_left_ATL_ant_tedana_D3(101:200,:),coords_acrossRun_left_ATL_ant_tedana_D3(201:300,:),coords_acrossRun_left_ATL_ant_tedana_D3(301:400,:)),3);
coords_acrossRun_left_ATL_pos_tedana_D1 = mean(cat(3,coords_acrossRun_left_ATL_pos_tedana_D1(1:100,:),coords_acrossRun_left_ATL_pos_tedana_D1(101:200,:),coords_acrossRun_left_ATL_pos_tedana_D1(201:300,:),coords_acrossRun_left_ATL_pos_tedana_D1(301:400,:)),3);
coords_acrossRun_left_ATL_pos_tedana_D2 = mean(cat(3,coords_acrossRun_left_ATL_pos_tedana_D2(1:100,:),coords_acrossRun_left_ATL_pos_tedana_D2(101:200,:),coords_acrossRun_left_ATL_pos_tedana_D2(201:300,:),coords_acrossRun_left_ATL_pos_tedana_D2(301:400,:)),3);
coords_acrossRun_left_ATL_pos_tedana_D3 = mean(cat(3,coords_acrossRun_left_ATL_pos_tedana_D3(1:100,:),coords_acrossRun_left_ATL_pos_tedana_D3(101:200,:),coords_acrossRun_left_ATL_pos_tedana_D3(201:300,:),coords_acrossRun_left_ATL_pos_tedana_D3(301:400,:)),3);
coords_acrossRun_right_ATL_ant_tedana_D1 = mean(cat(3,coords_acrossRun_right_ATL_ant_tedana_D1(1:100,:),coords_acrossRun_right_ATL_ant_tedana_D1(101:200,:),coords_acrossRun_right_ATL_ant_tedana_D1(201:300,:),coords_acrossRun_right_ATL_ant_tedana_D1(301:400,:)),3);
coords_acrossRun_right_ATL_ant_tedana_D2 = mean(cat(3,coords_acrossRun_right_ATL_ant_tedana_D2(1:100,:),coords_acrossRun_right_ATL_ant_tedana_D2(101:200,:),coords_acrossRun_right_ATL_ant_tedana_D2(201:300,:),coords_acrossRun_right_ATL_ant_tedana_D2(301:400,:)),3);
coords_acrossRun_right_ATL_ant_tedana_D3 = mean(cat(3,coords_acrossRun_right_ATL_ant_tedana_D3(1:100,:),coords_acrossRun_right_ATL_ant_tedana_D3(101:200,:),coords_acrossRun_right_ATL_ant_tedana_D3(201:300,:),coords_acrossRun_right_ATL_ant_tedana_D3(301:400,:)),3);
coords_acrossRun_right_ATL_pos_tedana_D1 = mean(cat(3,coords_acrossRun_right_ATL_pos_tedana_D1(1:100,:),coords_acrossRun_right_ATL_pos_tedana_D1(101:200,:),coords_acrossRun_right_ATL_pos_tedana_D1(201:300,:),coords_acrossRun_right_ATL_pos_tedana_D1(301:400,:)),3);
coords_acrossRun_right_ATL_pos_tedana_D2 = mean(cat(3,coords_acrossRun_right_ATL_pos_tedana_D2(1:100,:),coords_acrossRun_right_ATL_pos_tedana_D2(101:200,:),coords_acrossRun_right_ATL_pos_tedana_D2(201:300,:),coords_acrossRun_right_ATL_pos_tedana_D2(301:400,:)),3);
coords_acrossRun_right_ATL_pos_tedana_D3 = mean(cat(3,coords_acrossRun_right_ATL_pos_tedana_D3(1:100,:),coords_acrossRun_right_ATL_pos_tedana_D3(101:200,:),coords_acrossRun_right_ATL_pos_tedana_D3(201:300,:),coords_acrossRun_right_ATL_pos_tedana_D3(301:400,:)),3);

% next 20 figures show correlations with each dimension for all stimuli,
% just animate, and just inanimate; whole brain/ROI and subdivisions. No
% error bars for participants (because no folds), but we have error bars
% for group results (because multiple participants)

% Dimension 1, all stimuli, whole brain/whole ROI
figure;
% across run - whole brain - t2star
subplot(6,2,1)
bar(subcode,corr(coords_acrossRun_whole_brain_t2star_D1,U(:,1)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Across run - whole brain - t2star')
% across run - whole brain - tedana
subplot(6,2,3)
bar(subcode,corr(coords_acrossRun_whole_brain_tedana_D1,U(:,1)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - whole brain - tedana')
% across run - left ATL - t2star
subplot(6,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_t2star_D1,U(:,1)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Across run - left ATL - t2star')
% across run - left ATL - tedana
subplot(6,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_tedana_D1,U(:,1)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ATL - tedana')
% across run - right ATL - t2star
subplot(6,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_t2star_D1,U(:,1)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ATL - t2star')
% across run - right ATL - tedana
subplot(6,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_tedana_D1,U(:,1)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ATL - tedana')
% averaged - whole brain - t2star
subplot(6,2,2)
bar(subcode,corr(coords_averaged_whole_brain_t2star_D1,U(:,1)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Averaged - whole brain - t2star')
% averaged - whole brain - tedana
subplot(6,2,4)
bar(subcode,corr(coords_averaged_whole_brain_tedana_D1,U(:,1)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - whole brain - tedana')
% averaged - left ATL - t2star
subplot(6,2,6)
bar(subcode,corr(coords_averaged_left_ATL_t2star_D1,U(:,1)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Averaged - left ATL - t2star')
% averaged - left ATL - tedana
subplot(6,2,8)
bar(subcode,corr(coords_averaged_left_ATL_tedana_D1,U(:,1)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ATL - tedana')
% averaged - right ATL - t2star
subplot(6,2,10)
bar(subcode,corr(coords_averaged_right_ATL_t2star_D1,U(:,1)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ATL - t2star')
% averaged - right ATL - tedana
subplot(6,2,12)
bar(subcode,corr(coords_averaged_right_ATL_tedana_D1,U(:,1)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ATL - tedana')
sgtitle('Correlations with predicted coordinates - dimension 1')

% Dimension 2, all stimuli, whole brain/whole ROI
figure;
% across run - whole brain - t2star
subplot(6,2,1)
bar(subcode,corr(coords_acrossRun_whole_brain_t2star_D2,U(:,2)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Across run - whole brain - t2star')
% across run - whole brain - tedana
subplot(6,2,3)
bar(subcode,corr(coords_acrossRun_whole_brain_tedana_D2,U(:,2)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - whole brain - tedana')
% across run - left ATL - t2star
subplot(6,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_t2star_D2,U(:,2)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Across run - left ATL - t2star')
% across run - left ATL - tedana
subplot(6,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_tedana_D2,U(:,2)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ATL - tedana')
% across run - right ATL - t2star
subplot(6,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_t2star_D2,U(:,2)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ATL - t2star')
% across run - right ATL - tedana
subplot(6,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_tedana_D2,U(:,2)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ATL - tedana')
% averaged - whole brain - t2star
subplot(6,2,2)
bar(subcode,corr(coords_averaged_whole_brain_t2star_D2,U(:,2)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Averaged - whole brain - t2star')
% averaged - whole brain - tedana
subplot(6,2,4)
bar(subcode,corr(coords_averaged_whole_brain_tedana_D2,U(:,2)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - whole brain - tedana')
% averaged - left ATL - t2star
subplot(6,2,6)
bar(subcode,corr(coords_averaged_left_ATL_t2star_D2,U(:,2)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Averaged - left ATL - t2star')
% averaged - left ATL - tedana
subplot(6,2,8)
bar(subcode,corr(coords_averaged_left_ATL_tedana_D2,U(:,2)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ATL - tedana')
% averaged - right ATL - t2star
subplot(6,2,10)
bar(subcode,corr(coords_averaged_right_ATL_t2star_D2,U(:,2)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ATL - t2star')
% averaged - right ATL - tedana
subplot(6,2,12)
bar(subcode,corr(coords_averaged_right_ATL_tedana_D2,U(:,2)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ATL - tedana')
sgtitle('Correlations with predicted coordinates - dimension 2')

% Dimension 3, all stimuli, whole brain/whole ROI
figure;
% across run - whole brain - t2star
subplot(6,2,1)
bar(subcode,corr(coords_acrossRun_whole_brain_t2star_D3,U(:,3)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Across run - whole brain - t2star')
% across run - whole brain - tedana
subplot(6,2,3)
bar(subcode,corr(coords_acrossRun_whole_brain_tedana_D3,U(:,3)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - whole brain - tedana')
% across run - left ATL - t2star
subplot(6,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_t2star_D3,U(:,3)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Across run - left ATL - t2star')
% across run - left ATL - tedana
subplot(6,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_tedana_D3,U(:,3)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ATL - tedana')
% across run - right ATL - t2star
subplot(6,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_t2star_D3,U(:,3)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ATL - t2star')
% across run - right ATL - tedana
subplot(6,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_tedana_D3,U(:,3)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ATL - tedana')
% averaged - whole brain - t2star
subplot(6,2,2)
bar(subcode,corr(coords_averaged_whole_brain_t2star_D3,U(:,3)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Averaged - whole brain - t2star')
% averaged - whole brain - tedana
subplot(6,2,4)
bar(subcode,corr(coords_averaged_whole_brain_tedana_D3,U(:,3)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - whole brain - tedana')
% averaged - left ATL - t2star
subplot(6,2,6)
bar(subcode,corr(coords_averaged_left_ATL_t2star_D3,U(:,3)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Averaged - left ATL - t2star')
% averaged - left ATL - tedana
subplot(6,2,8)
bar(subcode,corr(coords_averaged_left_ATL_tedana_D3,U(:,3)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ATL - tedana')
% averaged - right ATL - t2star
subplot(6,2,10)
bar(subcode,corr(coords_averaged_right_ATL_t2star_D3,U(:,3)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ATL - t2star')
% averaged - right ATL - tedana
subplot(6,2,12)
bar(subcode,corr(coords_averaged_right_ATL_tedana_D3,U(:,3)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ATL - tedana')
sgtitle('Correlations with predicted coordinates - dimension 3')

% Dimension 1, animate only, whole brain/whole ROI
figure;
% across run - whole brain - t2star
subplot(6,2,1)
bar(subcode,corr(coords_acrossRun_whole_brain_t2star_D1(1:50,:),U(1:50,1)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Across run - whole brain - t2star')
% across run - whole brain - tedana
subplot(6,2,3)
bar(subcode,corr(coords_acrossRun_whole_brain_tedana_D1(1:50,:),U(1:50,1)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - whole brain - tedana')
% across run - left ATL - t2star
subplot(6,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_t2star_D1(1:50,:),U(1:50,1)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Across run - left ATL - t2star')
% across run - left ATL - tedana
subplot(6,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_tedana_D1(1:50,:),U(1:50,1)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ATL - tedana')
% across run - right ATL - t2star
subplot(6,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_t2star_D1(1:50,:),U(1:50,1)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ATL - t2star')
% across run - right ATL - tedana
subplot(6,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_tedana_D1(1:50,:),U(1:50,1)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ATL - tedana')
% averaged - whole brain - t2star
subplot(6,2,2)
bar(subcode,corr(coords_averaged_whole_brain_t2star_D1(1:50,:),U(1:50,1)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Averaged - whole brain - t2star')
% averaged - whole brain - tedana
subplot(6,2,4)
bar(subcode,corr(coords_averaged_whole_brain_tedana_D1(1:50,:),U(1:50,1)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - whole brain - tedana')
% averaged - left ATL - t2star
subplot(6,2,6)
bar(subcode,corr(coords_averaged_left_ATL_t2star_D1(1:50,:),U(1:50,1)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Averaged - left ATL - t2star')
% averaged - left ATL - tedana
subplot(6,2,8)
bar(subcode,corr(coords_averaged_left_ATL_tedana_D1(1:50,:),U(1:50,1)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ATL - tedana')
% averaged - right ATL - t2star
subplot(6,2,10)
bar(subcode,corr(coords_averaged_right_ATL_t2star_D1(1:50,:),U(1:50,1)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ATL - t2star')
% averaged - right ATL - tedana
subplot(6,2,12)
bar(subcode,corr(coords_averaged_right_ATL_tedana_D1(1:50,:),U(1:50,1)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ATL - tedana')
sgtitle('Correlations with predicted coordinates - animate only - dimension 1')

% Dimension 2, animate only, whole brain/whole ROI
figure;
% across run - whole brain - t2star
subplot(6,2,1)
bar(subcode,corr(coords_acrossRun_whole_brain_t2star_D2(1:50,:),U(1:50,2)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Across run - whole brain - t2star')
% across run - whole brain - tedana
subplot(6,2,3)
bar(subcode,corr(coords_acrossRun_whole_brain_tedana_D2(1:50,:),U(1:50,2)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - whole brain - tedana')
% across run - left ATL - t2star
subplot(6,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_t2star_D2(1:50,:),U(1:50,2)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Across run - left ATL - t2star')
% across run - left ATL - tedana
subplot(6,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_tedana_D2(1:50,:),U(1:50,2)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ATL - tedana')
% across run - right ATL - t2star
subplot(6,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_t2star_D2(1:50,:),U(1:50,2)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ATL - t2star')
% across run - right ATL - tedana
subplot(6,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_tedana_D2(1:50,:),U(1:50,2)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ATL - tedana')
% averaged - whole brain - t2star
subplot(6,2,2)
bar(subcode,corr(coords_averaged_whole_brain_t2star_D2(1:50,:),U(1:50,2)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Averaged - whole brain - t2star')
% averaged - whole brain - tedana
subplot(6,2,4)
bar(subcode,corr(coords_averaged_whole_brain_tedana_D2(1:50,:),U(1:50,2)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - whole brain - tedana')
% averaged - left ATL - t2star
subplot(6,2,6)
bar(subcode,corr(coords_averaged_left_ATL_t2star_D2(1:50,:),U(1:50,2)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Averaged - left ATL - t2star')
% averaged - left ATL - tedana
subplot(6,2,8)
bar(subcode,corr(coords_averaged_left_ATL_tedana_D2(1:50,:),U(1:50,2)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ATL - tedana')
% averaged - right ATL - t2star
subplot(6,2,10)
bar(subcode,corr(coords_averaged_right_ATL_t2star_D2(1:50,:),U(1:50,2)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ATL - t2star')
% averaged - right ATL - tedana
subplot(6,2,12)
bar(subcode,corr(coords_averaged_right_ATL_tedana_D2(1:50,:),U(1:50,2)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ATL - tedana')
sgtitle('Correlations with predicted coordinates - animate only - dimension 2')

% Dimension 3, animate only, whole brain/whole ROI
figure;
% across run - whole brain - t2star
subplot(6,2,1)
bar(subcode,corr(coords_acrossRun_whole_brain_t2star_D3(1:50,:),U(1:50,3)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Across run - whole brain - t2star')
% across run - whole brain - tedana
subplot(6,2,3)
bar(subcode,corr(coords_acrossRun_whole_brain_tedana_D3(1:50,:),U(1:50,3)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - whole brain - tedana')
% across run - left ATL - t2star
subplot(6,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_t2star_D3(1:50,:),U(1:50,3)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Across run - left ATL - t2star')
% across run - left ATL - tedana
subplot(6,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_tedana_D3(1:50,:),U(1:50,3)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ATL - tedana')
% across run - right ATL - t2star
subplot(6,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_t2star_D3(1:50,:),U(1:50,3)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ATL - t2star')
% across run - right ATL - tedana
subplot(6,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_tedana_D3(1:50,:),U(1:50,3)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ATL - tedana')
% averaged - whole brain - t2star
subplot(6,2,2)
bar(subcode,corr(coords_averaged_whole_brain_t2star_D3(1:50,:),U(1:50,3)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Averaged - whole brain - t2star')
% averaged - whole brain - tedana
subplot(6,2,4)
bar(subcode,corr(coords_averaged_whole_brain_tedana_D3(1:50,:),U(1:50,3)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - whole brain - tedana')
% averaged - left ATL - t2star
subplot(6,2,6)
bar(subcode,corr(coords_averaged_left_ATL_t2star_D3(1:50,:),U(1:50,3)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Averaged - left ATL - t2star')
% averaged - left ATL - tedana
subplot(6,2,8)
bar(subcode,corr(coords_averaged_left_ATL_tedana_D3(1:50,:),U(1:50,3)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ATL - tedana')
% averaged - right ATL - t2star
subplot(6,2,10)
bar(subcode,corr(coords_averaged_right_ATL_t2star_D3(1:50,:),U(1:50,3)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ATL - t2star')
% averaged - right ATL - tedana
subplot(6,2,12)
bar(subcode,corr(coords_averaged_right_ATL_tedana_D3(1:50,:),U(1:50,3)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ATL - tedana')
sgtitle('Correlations with predicted coordinates - animate only - dimension 3')

% Dimension 1, inanimate only, whole brain/whole ROI
figure;
% across run - whole brain - t2star
subplot(6,2,1)
bar(subcode,corr(coords_acrossRun_whole_brain_t2star_D1(51:100,:),U(51:100,1)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Across run - whole brain - t2star')
% across run - whole brain - tedana
subplot(6,2,3)
bar(subcode,corr(coords_acrossRun_whole_brain_tedana_D1(51:100,:),U(51:100,1)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - whole brain - tedana')
% across run - left ATL - t2star
subplot(6,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_t2star_D1(51:100,:),U(51:100,1)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Across run - left ATL - t2star')
% across run - left ATL - tedana
subplot(6,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_tedana_D1(51:100,:),U(51:100,1)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ATL - tedana')
% across run - right ATL - t2star
subplot(6,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_t2star_D1(51:100,:),U(51:100,1)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ATL - t2star')
% across run - right ATL - tedana
subplot(6,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_tedana_D1(51:100,:),U(51:100,1)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ATL - tedana')
% averaged - whole brain - t2star
subplot(6,2,2)
bar(subcode,corr(coords_averaged_whole_brain_t2star_D1(51:100,:),U(51:100,1)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Averaged - whole brain - t2star')
% averaged - whole brain - tedana
subplot(6,2,4)
bar(subcode,corr(coords_averaged_whole_brain_tedana_D1(51:100,:),U(51:100,1)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - whole brain - tedana')
% averaged - left ATL - t2star
subplot(6,2,6)
bar(subcode,corr(coords_averaged_left_ATL_t2star_D1(51:100,:),U(51:100,1)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Averaged - left ATL - t2star')
% averaged - left ATL - tedana
subplot(6,2,8)
bar(subcode,corr(coords_averaged_left_ATL_tedana_D1(51:100,:),U(51:100,1)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ATL - tedana')
% averaged - right ATL - t2star
subplot(6,2,10)
bar(subcode,corr(coords_averaged_right_ATL_t2star_D1(51:100,:),U(51:100,1)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ATL - t2star')
% averaged - right ATL - tedana
subplot(6,2,12)
bar(subcode,corr(coords_averaged_right_ATL_tedana_D1(51:100,:),U(51:100,1)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ATL - tedana')
sgtitle('Correlations with predicted coordinates - inanimate only - dimension 1')

% Dimension 2, inanimate only, whole brain/whole ROI
figure;
% across run - whole brain - t2star
subplot(6,2,1)
bar(subcode,corr(coords_acrossRun_whole_brain_t2star_D2(51:100,:),U(51:100,2)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Across run - whole brain - t2star')
% across run - whole brain - tedana
subplot(6,2,3)
bar(subcode,corr(coords_acrossRun_whole_brain_tedana_D2(51:100,:),U(51:100,2)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - whole brain - tedana')
% across run - left ATL - t2star
subplot(6,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_t2star_D2(51:100,:),U(51:100,2)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Across run - left ATL - t2star')
% across run - left ATL - tedana
subplot(6,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_tedana_D2(51:100,:),U(51:100,2)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ATL - tedana')
% across run - right ATL - t2star
subplot(6,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_t2star_D2(51:100,:),U(51:100,2)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ATL - t2star')
% across run - right ATL - tedana
subplot(6,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_tedana_D2(51:100,:),U(51:100,2)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ATL - tedana')
% averaged - whole brain - t2star
subplot(6,2,2)
bar(subcode,corr(coords_averaged_whole_brain_t2star_D2(51:100,:),U(51:100,2)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Averaged - whole brain - t2star')
% averaged - whole brain - tedana
subplot(6,2,4)
bar(subcode,corr(coords_averaged_whole_brain_tedana_D2(51:100,:),U(51:100,2)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - whole brain - tedana')
% averaged - left ATL - t2star
subplot(6,2,6)
bar(subcode,corr(coords_averaged_left_ATL_t2star_D2(51:100,:),U(51:100,2)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Averaged - left ATL - t2star')
% averaged - left ATL - tedana
subplot(6,2,8)
bar(subcode,corr(coords_averaged_left_ATL_tedana_D2(51:100,:),U(51:100,2)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ATL - tedana')
% averaged - right ATL - t2star
subplot(6,2,10)
bar(subcode,corr(coords_averaged_right_ATL_t2star_D2(51:100,:),U(51:100,2)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ATL - t2star')
% averaged - right ATL - tedana
subplot(6,2,12)
bar(subcode,corr(coords_averaged_right_ATL_tedana_D2(51:100,:),U(51:100,2)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ATL - tedana')
sgtitle('Correlations with predicted coordinates - inanimate only - dimension 2')

% Dimension 3, inanimate only, whole brain/whole ROI
figure;
% across run - whole brain - t2star
subplot(6,2,1)
bar(subcode,corr(coords_acrossRun_whole_brain_t2star_D3(51:100,:),U(51:100,3)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Across run - whole brain - t2star')
% across run - whole brain - tedana
subplot(6,2,3)
bar(subcode,corr(coords_acrossRun_whole_brain_tedana_D3(51:100,:),U(51:100,3)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - whole brain - tedana')
% across run - left ATL - t2star
subplot(6,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_t2star_D3(51:100,:),U(51:100,3)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Across run - left ATL - t2star')
% across run - left ATL - tedana
subplot(6,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_tedana_D3(51:100,:),U(51:100,3)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ATL - tedana')
% across run - right ATL - t2star
subplot(6,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_t2star_D3(51:100,:),U(51:100,3)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ATL - t2star')
% across run - right ATL - tedana
subplot(6,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_tedana_D3(51:100,:),U(51:100,3)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ATL - tedana')
% averaged - whole brain - t2star
subplot(6,2,2)
bar(subcode,corr(coords_averaged_whole_brain_t2star_D3(51:100,:),U(51:100,3)),'FaceColor','#77AC30')
ylim([-0.1,1])
title('Averaged - whole brain - t2star')
% averaged - whole brain - tedana
subplot(6,2,4)
bar(subcode,corr(coords_averaged_whole_brain_tedana_D3(51:100,:),U(51:100,3)),'FaceColor','#77AC30','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - whole brain - tedana')
% averaged - left ATL - t2star
subplot(6,2,6)
bar(subcode,corr(coords_averaged_left_ATL_t2star_D3(51:100,:),U(51:100,3)),'FaceColor','#A2142F')
ylim([-0.1,1])
title('Averaged - left ATL - t2star')
% averaged - left ATL - tedana
subplot(6,2,8)
bar(subcode,corr(coords_averaged_left_ATL_tedana_D3(51:100,:),U(51:100,3)),'FaceColor','#A2142F','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ATL - tedana')
% averaged - right ATL - t2star
subplot(6,2,10)
bar(subcode,corr(coords_averaged_right_ATL_t2star_D3(51:100,:),U(51:100,3)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ATL - t2star')
% averaged - right ATL - tedana
subplot(6,2,12)
bar(subcode,corr(coords_averaged_right_ATL_tedana_D3(51:100,:),U(51:100,3)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ATL - tedana')
sgtitle('Correlations with predicted coordinates - inanimate only - dimension 3')

% group results, all stimuli, whole brain/whole ROI
groupcorrall = zeros(size(subcode,2),36);
groupcorrall(:,1) = corr(coords_acrossRun_whole_brain_t2star_D1,U(:,1));
groupcorrall(:,2) = corr(coords_acrossRun_whole_brain_t2star_D2,U(:,2));
groupcorrall(:,3) = corr(coords_acrossRun_whole_brain_t2star_D3,U(:,3));
groupcorrall(:,4) = corr(coords_acrossRun_whole_brain_tedana_D1,U(:,1));
groupcorrall(:,5) = corr(coords_acrossRun_whole_brain_tedana_D2,U(:,2));
groupcorrall(:,6) = corr(coords_acrossRun_whole_brain_tedana_D3,U(:,3));
groupcorrall(:,7) = corr(coords_acrossRun_left_ATL_t2star_D1,U(:,1));
groupcorrall(:,8) = corr(coords_acrossRun_left_ATL_t2star_D2,U(:,2));
groupcorrall(:,9) = corr(coords_acrossRun_left_ATL_t2star_D3,U(:,3));
groupcorrall(:,10) = corr(coords_acrossRun_left_ATL_tedana_D1,U(:,1));
groupcorrall(:,11) = corr(coords_acrossRun_left_ATL_tedana_D2,U(:,2));
groupcorrall(:,12) = corr(coords_acrossRun_left_ATL_tedana_D3,U(:,3));
groupcorrall(:,13) = corr(coords_acrossRun_right_ATL_t2star_D1,U(:,1));
groupcorrall(:,14) = corr(coords_acrossRun_right_ATL_t2star_D2,U(:,2));
groupcorrall(:,15) = corr(coords_acrossRun_right_ATL_t2star_D3,U(:,3));
groupcorrall(:,16) = corr(coords_acrossRun_right_ATL_tedana_D1,U(:,1));
groupcorrall(:,17) = corr(coords_acrossRun_right_ATL_tedana_D2,U(:,2));
groupcorrall(:,18) = corr(coords_averaged_right_ATL_tedana_D3,U(:,3));
groupcorrall(:,19) = corr(coords_averaged_whole_brain_t2star_D1,U(:,1));
groupcorrall(:,20) = corr(coords_averaged_whole_brain_t2star_D2,U(:,2));
groupcorrall(:,21) = corr(coords_averaged_whole_brain_t2star_D3,U(:,3));
groupcorrall(:,22) = corr(coords_averaged_whole_brain_tedana_D1,U(:,1));
groupcorrall(:,23) = corr(coords_averaged_whole_brain_tedana_D2,U(:,2));
groupcorrall(:,24) = corr(coords_averaged_whole_brain_tedana_D3,U(:,3));
groupcorrall(:,25) = corr(coords_averaged_left_ATL_t2star_D1,U(:,1));
groupcorrall(:,26) = corr(coords_averaged_left_ATL_t2star_D2,U(:,2));
groupcorrall(:,27) = corr(coords_averaged_left_ATL_t2star_D3,U(:,3));
groupcorrall(:,28) = corr(coords_averaged_left_ATL_tedana_D1,U(:,1));
groupcorrall(:,29) = corr(coords_averaged_left_ATL_tedana_D2,U(:,2));
groupcorrall(:,30) = corr(coords_averaged_left_ATL_tedana_D3,U(:,3));
groupcorrall(:,31) = corr(coords_averaged_right_ATL_t2star_D1,U(:,1));
groupcorrall(:,32) = corr(coords_averaged_right_ATL_t2star_D2,U(:,2));
groupcorrall(:,33) = corr(coords_averaged_right_ATL_t2star_D3,U(:,3));
groupcorrall(:,34) = corr(coords_averaged_right_ATL_tedana_D1,U(:,1));
groupcorrall(:,35) = corr(coords_averaged_right_ATL_tedana_D2,U(:,2));
groupcorrall(:,36) = corr(coords_averaged_right_ATL_tedana_D3,U(:,3));
colnames = {'AR-WB-t2star-D1','AR-WB-t2star-D2','AR-WB-t2star-D3','AR-WB-tedana-D1','AR-WB-tedana-D2','AR-WB-tedana-D3','AR-lATL-t2star-D1','AR-lATL-t2star-D2','AR-lATL-t2star-D3','AR-lATL-tedana-D1','AR-lATL-tedana-D2','AR-lATL-tedana-D3','AR-rATL-t2star-D1','AR-rATL-t2star-D2','AR-rATL-t2star-D3','AR-rATL-tedana-D1','AR-rATL-tedana-D2','AR-rATL-tedana-D3','av-WB-t2star-D1','av-WB-t2star-D2','av-WB-t2star-D3','av-WB-tedana-D1','av-WB-tedana-D2','av-WB-tedana-D3','av-lATL-t2star-D1','av-lATL-t2star-D2','av-lATL-t2star-D3','av-lATL-tedana-D1','av-lATL-tedana-D2','av-lATL-tedana-D3','av-rATL-t2star-D1','av-rATL-t2star-D2','av-rATL-t2star-D3','av-rATL-tedana-D1','av-rATL-tedana-D2','av-rATL-tedana-D3'};

figure;
b = bar(colnames,mean(groupcorrall),'FaceColor','flat');
% normalised RGB triplets
colours = [119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189;119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189]/255;
b.CData = colours;
ylim([-0.1,1])
title('Group results - correlations with predicted coordinates')
hold on
confint = 1.96*(std(groupcorrall)/sqrt(size(groupcorrall,1)));
errorbar(1:size(groupcorrall,2),mean(groupcorrall,1),confint,'LineStyle','none','Color','black');

% group results, animate only, whole brain/whole ROI
groupcorranimate = zeros(size(subcode,2),36);
groupcorranimate(:,1) = corr(coords_acrossRun_whole_brain_t2star_D1(1:50,:),U(1:50,1));
groupcorranimate(:,2) = corr(coords_acrossRun_whole_brain_t2star_D2(1:50,:),U(1:50,2));
groupcorranimate(:,3) = corr(coords_acrossRun_whole_brain_t2star_D3(1:50,:),U(1:50,3));
groupcorranimate(:,4) = corr(coords_acrossRun_whole_brain_tedana_D1(1:50,:),U(1:50,1));
groupcorranimate(:,5) = corr(coords_acrossRun_whole_brain_tedana_D2(1:50,:),U(1:50,2));
groupcorranimate(:,6) = corr(coords_acrossRun_whole_brain_tedana_D3(1:50,:),U(1:50,3));
groupcorranimate(:,7) = corr(coords_acrossRun_left_ATL_t2star_D1(1:50,:),U(1:50,1));
groupcorranimate(:,8) = corr(coords_acrossRun_left_ATL_t2star_D2(1:50,:),U(1:50,2));
groupcorranimate(:,9) = corr(coords_acrossRun_left_ATL_t2star_D3(1:50,:),U(1:50,3));
groupcorranimate(:,10) = corr(coords_acrossRun_left_ATL_tedana_D1(1:50,:),U(1:50,1));
groupcorranimate(:,11) = corr(coords_acrossRun_left_ATL_tedana_D2(1:50,:),U(1:50,2));
groupcorranimate(:,12) = corr(coords_acrossRun_left_ATL_tedana_D3(1:50,:),U(1:50,3));
groupcorranimate(:,13) = corr(coords_acrossRun_right_ATL_t2star_D1(1:50,:),U(1:50,1));
groupcorranimate(:,14) = corr(coords_acrossRun_right_ATL_t2star_D2(1:50,:),U(1:50,2));
groupcorranimate(:,15) = corr(coords_acrossRun_right_ATL_t2star_D3(1:50,:),U(1:50,3));
groupcorranimate(:,16) = corr(coords_acrossRun_right_ATL_tedana_D1(1:50,:),U(1:50,1));
groupcorranimate(:,17) = corr(coords_acrossRun_right_ATL_tedana_D2(1:50,:),U(1:50,2));
groupcorranimate(:,18) = corr(coords_acrossRun_right_ATL_tedana_D3(1:50,:),U(1:50,3));
groupcorranimate(:,19) = corr(coords_averaged_whole_brain_t2star_D1(1:50,:),U(1:50,1));
groupcorranimate(:,20) = corr(coords_averaged_whole_brain_t2star_D2(1:50,:),U(1:50,2));
groupcorranimate(:,21) = corr(coords_averaged_whole_brain_t2star_D3(1:50,:),U(1:50,3));
groupcorranimate(:,22) = corr(coords_averaged_whole_brain_tedana_D1(1:50,:),U(1:50,1));
groupcorranimate(:,23) = corr(coords_averaged_whole_brain_tedana_D2(1:50,:),U(1:50,2));
groupcorranimate(:,24) = corr(coords_averaged_whole_brain_tedana_D3(1:50,:),U(1:50,3));
groupcorranimate(:,25) = corr(coords_averaged_left_ATL_t2star_D1(1:50,:),U(1:50,1));
groupcorranimate(:,26) = corr(coords_averaged_left_ATL_t2star_D2(1:50,:),U(1:50,2));
groupcorranimate(:,27) = corr(coords_averaged_left_ATL_t2star_D3(1:50,:),U(1:50,3));
groupcorranimate(:,28) = corr(coords_averaged_left_ATL_tedana_D1(1:50,:),U(1:50,1));
groupcorranimate(:,29) = corr(coords_averaged_left_ATL_tedana_D2(1:50,:),U(1:50,2));
groupcorranimate(:,30) = corr(coords_averaged_left_ATL_tedana_D3(1:50,:),U(1:50,3));
groupcorranimate(:,31) = corr(coords_averaged_right_ATL_t2star_D1(1:50,:),U(1:50,1));
groupcorranimate(:,32) = corr(coords_averaged_right_ATL_t2star_D2(1:50,:),U(1:50,2));
groupcorranimate(:,33) = corr(coords_averaged_right_ATL_t2star_D3(1:50,:),U(1:50,3));
groupcorranimate(:,34) = corr(coords_averaged_right_ATL_tedana_D1(1:50,:),U(1:50,1));
groupcorranimate(:,35) = corr(coords_averaged_right_ATL_tedana_D2(1:50,:),U(1:50,2));
groupcorranimate(:,36) = corr(coords_averaged_right_ATL_tedana_D3(1:50,:),U(1:50,3));
colnames = {'AR-WB-t2star-D1','AR-WB-t2star-D2','AR-WB-t2star-D3','AR-WB-tedana-D1','AR-WB-tedana-D2','AR-WB-tedana-D3','AR-lATL-t2star-D1','AR-lATL-t2star-D2','AR-lATL-t2star-D3','AR-lATL-tedana-D1','AR-lATL-tedana-D2','AR-lATL-tedana-D3','AR-rATL-t2star-D1','AR-rATL-t2star-D2','AR-rATL-t2star-D3','AR-rATL-tedana-D1','AR-rATL-tedana-D2','AR-rATL-tedana-D3','av-WB-t2star-D1','av-WB-t2star-D2','av-WB-t2star-D3','av-WB-tedana-D1','av-WB-tedana-D2','av-WB-tedana-D3','av-lATL-t2star-D1','av-lATL-t2star-D2','av-lATL-t2star-D3','av-lATL-tedana-D1','av-lATL-tedana-D2','av-lATL-tedana-D3','av-rATL-t2star-D1','av-rATL-t2star-D2','av-rATL-t2star-D3','av-rATL-tedana-D1','av-rATL-tedana-D2','av-rATL-tedana-D3'};

figure;
b = bar(colnames,mean(groupcorranimate),'FaceColor','flat');
% normalised RGB triplets
colours = [119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189;119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189]/255;
b.CData = colours;
ylim([-0.1,1])
title('Group results - animate only - correlations with predicted coordinates')
hold on
confint = 1.96*(std(groupcorranimate)/sqrt(size(groupcorranimate,1)));
errorbar(1:size(groupcorranimate,2),mean(groupcorranimate,1),confint,'LineStyle','none','Color','black');

% group results, inanimate only, whole brain/whole ROI
groupcorrinanimate = zeros(size(subcode,2),36);
groupcorrinanimate(:,1) = corr(coords_acrossRun_whole_brain_t2star_D1(51:100,:),U(51:100,1));
groupcorrinanimate(:,2) = corr(coords_acrossRun_whole_brain_t2star_D2(51:100,:),U(51:100,2));
groupcorrinanimate(:,3) = corr(coords_acrossRun_whole_brain_t2star_D3(51:100,:),U(51:100,3));
groupcorrinanimate(:,4) = corr(coords_acrossRun_whole_brain_tedana_D1(51:100,:),U(51:100,1));
groupcorrinanimate(:,5) = corr(coords_acrossRun_whole_brain_tedana_D2(51:100,:),U(51:100,2));
groupcorrinanimate(:,6) = corr(coords_acrossRun_whole_brain_tedana_D3(51:100,:),U(51:100,3));
groupcorrinanimate(:,7) = corr(coords_acrossRun_left_ATL_t2star_D1(51:100,:),U(51:100,1));
groupcorrinanimate(:,8) = corr(coords_acrossRun_left_ATL_t2star_D2(51:100,:),U(51:100,2));
groupcorrinanimate(:,9) = corr(coords_acrossRun_left_ATL_t2star_D3(51:100,:),U(51:100,3));
groupcorrinanimate(:,10) = corr(coords_acrossRun_left_ATL_tedana_D1(51:100,:),U(51:100,1));
groupcorrinanimate(:,11) = corr(coords_acrossRun_left_ATL_tedana_D2(51:100,:),U(51:100,2));
groupcorrinanimate(:,12) = corr(coords_acrossRun_left_ATL_tedana_D3(51:100,:),U(51:100,3));
groupcorrinanimate(:,13) = corr(coords_acrossRun_right_ATL_t2star_D1(51:100,:),U(51:100,1));
groupcorrinanimate(:,14) = corr(coords_acrossRun_right_ATL_t2star_D2(51:100,:),U(51:100,2));
groupcorrinanimate(:,15) = corr(coords_acrossRun_right_ATL_t2star_D3(51:100,:),U(51:100,3));
groupcorrinanimate(:,16) = corr(coords_acrossRun_right_ATL_tedana_D1(51:100,:),U(51:100,1));
groupcorrinanimate(:,17) = corr(coords_acrossRun_right_ATL_tedana_D2(51:100,:),U(51:100,2));
groupcorrinanimate(:,18) = corr(coords_acrossRun_right_ATL_tedana_D3(51:100,:),U(51:100,3));
groupcorrinanimate(:,19) = corr(coords_averaged_whole_brain_t2star_D1(51:100,:),U(51:100,1));
groupcorrinanimate(:,20) = corr(coords_averaged_whole_brain_t2star_D2(51:100,:),U(51:100,2));
groupcorrinanimate(:,21) = corr(coords_averaged_whole_brain_t2star_D3(51:100,:),U(51:100,3));
groupcorrinanimate(:,22) = corr(coords_averaged_whole_brain_tedana_D1(51:100,:),U(51:100,1));
groupcorrinanimate(:,23) = corr(coords_averaged_whole_brain_tedana_D2(51:100,:),U(51:100,2));
groupcorrinanimate(:,24) = corr(coords_averaged_whole_brain_tedana_D3(51:100,:),U(51:100,3));
groupcorrinanimate(:,25) = corr(coords_averaged_left_ATL_t2star_D1(51:100,:),U(51:100,1));
groupcorrinanimate(:,26) = corr(coords_averaged_left_ATL_t2star_D2(51:100,:),U(51:100,2));
groupcorrinanimate(:,27) = corr(coords_averaged_left_ATL_t2star_D3(51:100,:),U(51:100,3));
groupcorrinanimate(:,28) = corr(coords_averaged_left_ATL_tedana_D1(51:100,:),U(51:100,1));
groupcorrinanimate(:,29) = corr(coords_averaged_left_ATL_tedana_D2(51:100,:),U(51:100,2));
groupcorrinanimate(:,30) = corr(coords_averaged_left_ATL_tedana_D3(51:100,:),U(51:100,3));
groupcorrinanimate(:,31) = corr(coords_averaged_right_ATL_t2star_D1(51:100,:),U(51:100,1));
groupcorrinanimate(:,32) = corr(coords_averaged_right_ATL_t2star_D2(51:100,:),U(51:100,2));
groupcorrinanimate(:,33) = corr(coords_averaged_right_ATL_t2star_D3(51:100,:),U(51:100,3));
groupcorrinanimate(:,34) = corr(coords_averaged_right_ATL_tedana_D1(51:100,:),U(51:100,1));
groupcorrinanimate(:,35) = corr(coords_averaged_right_ATL_tedana_D2(51:100,:),U(51:100,2));
groupcorrinanimate(:,36) = corr(coords_averaged_right_ATL_tedana_D3(51:100,:),U(51:100,3));
colnames = {'AR-WB-t2star-D1','AR-WB-t2star-D2','AR-WB-t2star-D3','AR-WB-tedana-D1','AR-WB-tedana-D2','AR-WB-tedana-D3','AR-lATL-t2star-D1','AR-lATL-t2star-D2','AR-lATL-t2star-D3','AR-lATL-tedana-D1','AR-lATL-tedana-D2','AR-lATL-tedana-D3','AR-rATL-t2star-D1','AR-rATL-t2star-D2','AR-rATL-t2star-D3','AR-rATL-tedana-D1','AR-rATL-tedana-D2','AR-rATL-tedana-D3','av-WB-t2star-D1','av-WB-t2star-D2','av-WB-t2star-D3','av-WB-tedana-D1','av-WB-tedana-D2','av-WB-tedana-D3','av-lATL-t2star-D1','av-lATL-t2star-D2','av-lATL-t2star-D3','av-lATL-tedana-D1','av-lATL-tedana-D2','av-lATL-tedana-D3','av-rATL-t2star-D1','av-rATL-t2star-D2','av-rATL-t2star-D3','av-rATL-tedana-D1','av-rATL-tedana-D2','av-rATL-tedana-D3'};

figure;
b = bar(colnames,mean(groupcorrinanimate),'FaceColor','flat');
% normalised RGB triplets
colours = [119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189;119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;119,172,48;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;162,20,47;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189;0,114,189]/255;
b.CData = colours;
ylim([-0.1,1])
title('Group results - inanimate only - correlations with predicted coordinates')
hold on
confint = 1.96*(std(groupcorrinanimate)/sqrt(size(groupcorrinanimate,1)));
errorbar(1:size(groupcorrinanimate,2),mean(groupcorrinanimate,1),confint,'LineStyle','none','Color','black');

% Dimension 1, all stimuli, anterior/posterior halves of ROI
figure;
% across run - left ROI (anterior) - t2star
subplot(8,2,1)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_t2star_D1,U(:,1)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Across run - left ROI (anterior) - t2star')
% across run - left ROI (anterior) - tedana
subplot(8,2,3)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_tedana_D1,U(:,1)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - left ROI (posterior) - t2star
subplot(8,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_t2star_D1,U(:,1)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Across run - left ROI (posterior) - t2star')
% across run - left ROI (posterior) - tedana
subplot(8,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_tedana_D1,U(:,1)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - right ROI (anterior) - t2star
subplot(8,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_t2star_D1,U(:,1)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ROI (anterior) - t2star')
% across run - right ROI (anterior) - tedana
subplot(8,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_tedana_D1,U(:,1)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% across run - right ROI (posterior) - t2star
subplot(8,2,13)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_t2star_D1,U(:,1)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Across run - right ROI (posterior) - t2star')
% across run - right ROI (posterior) - tedana
subplot(8,2,15)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_tedana_D1,U(:,1)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% averaged - left ROI (anterior) - t2star
subplot(8,2,2)
bar(subcode,corr(coords_averaged_left_ATL_ant_t2star_D1,U(:,1)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - t2star')
% averaged - left ROI (anterior) - tedana
subplot(8,2,4)
bar(subcode,corr(coords_averaged_left_ATL_ant_tedana_D1,U(:,1)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - tedana')
% averaged - left ROI (posterior) - t2star
subplot(8,2,6)
bar(subcode,corr(coords_averaged_left_ATL_pos_t2star_D1,U(:,1)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - t2star')
% averaged - left ROI (posterior) - tedana
subplot(8,2,8)
bar(subcode,corr(coords_averaged_left_ATL_pos_tedana_D1,U(:,1)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - tedana')
% averaged - right ROI (anterior) - t2star
subplot(8,2,10)
bar(subcode,corr(coords_averaged_right_ATL_ant_t2star_D1,U(:,1)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - t2star')
% averaged - right ROI (anterior) - tedana
subplot(8,2,12)
bar(subcode,corr(coords_averaged_right_ATL_ant_tedana_D1,U(:,1)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
hold on
% averaged - right ROI (posterior) - t2star
subplot(8,2,14)
bar(subcode,corr(coords_averaged_right_ATL_pos_t2star_D1,U(:,1)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Averaged - right ROI (posterior) - t2star')
% averaged - right ROI (posterior) - tedana
subplot(8,2,16)
bar(subcode,corr(coords_averaged_right_ATL_pos_tedana_D1,U(:,1)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
sgtitle('Correlations with predicted coordinates - all stimuli - dimension 1 - anterior and posterior halves of ROI')

% Dimension 2, all stimuli, anterior/posterior halves of ROI
figure;
% across run - left ROI (anterior) - t2star
subplot(8,2,1)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_t2star_D2,U(:,2)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Across run - left ROI (anterior) - t2star')
% across run - left ROI (anterior) - tedana
subplot(8,2,3)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_tedana_D2,U(:,2)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - left ROI (posterior) - t2star
subplot(8,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_t2star_D2,U(:,2)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Across run - left ROI (posterior) - t2star')
% across run - left ROI (posterior) - tedana
subplot(8,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_tedana_D2,U(:,2)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - right ROI (anterior) - t2star
subplot(8,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_t2star_D2,U(:,2)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ROI (anterior) - t2star')
% across run - right ROI (anterior) - tedana
subplot(8,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_tedana_D2,U(:,2)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% across run - right ROI (posterior) - t2star
subplot(8,2,13)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_t2star_D2,U(:,2)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Across run - right ROI (posterior) - t2star')
% across run - right ROI (posterior) - tedana
subplot(8,2,15)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_tedana_D2,U(:,2)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% averaged - left ROI (anterior) - t2star
subplot(8,2,2)
bar(subcode,corr(coords_averaged_left_ATL_ant_t2star_D2,U(:,2)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - t2star')
% averaged - left ROI (anterior) - tedana
subplot(8,2,4)
bar(subcode,corr(coords_averaged_left_ATL_ant_tedana_D2,U(:,2)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - tedana')
% averaged - left ROI (posterior) - t2star
subplot(8,2,6)
bar(subcode,corr(coords_averaged_left_ATL_pos_t2star_D2,U(:,2)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - t2star')
% averaged - left ROI (posterior) - tedana
subplot(8,2,8)
bar(subcode,corr(coords_averaged_left_ATL_pos_tedana_D2,U(:,2)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - tedana')
% averaged - right ROI (anterior) - t2star
subplot(8,2,10)
bar(subcode,corr(coords_averaged_right_ATL_ant_t2star_D2,U(:,2)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - t2star')
% averaged - right ROI (anterior) - tedana
subplot(8,2,12)
bar(subcode,corr(coords_averaged_right_ATL_ant_tedana_D2,U(:,2)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
hold on
% averaged - right ROI (posterior) - t2star
subplot(8,2,14)
bar(subcode,corr(coords_averaged_right_ATL_pos_t2star_D2,U(:,2)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Averaged - right ROI (posterior) - t2star')
% averaged - right ROI (posterior) - tedana
subplot(8,2,16)
bar(subcode,corr(coords_averaged_right_ATL_pos_tedana_D2,U(:,2)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
sgtitle('Correlations with predicted coordinates - all stimuli - dimension 2 - anterior and posterior halves of ROI')

% Dimension 3, all stimuli, anterior/posterior halves of ROI
figure;
% across run - left ROI (anterior) - t2star
subplot(8,2,1)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_t2star_D3,U(:,3)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Across run - left ROI (anterior) - t2star')
% across run - left ROI (anterior) - tedana
subplot(8,2,3)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_tedana_D3,U(:,3)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - left ROI (posterior) - t2star
subplot(8,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_t2star_D3,U(:,3)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Across run - left ROI (posterior) - t2star')
% across run - left ROI (posterior) - tedana
subplot(8,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_tedana_D3,U(:,3)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - right ROI (anterior) - t2star
subplot(8,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_t2star_D3,U(:,3)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ROI (anterior) - t2star')
% across run - right ROI (anterior) - tedana
subplot(8,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_tedana_D3,U(:,3)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% across run - right ROI (posterior) - t2star
subplot(8,2,13)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_t2star_D3,U(:,3)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Across run - right ROI (posterior) - t2star')
% across run - right ROI (posterior) - tedana
subplot(8,2,15)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_tedana_D3,U(:,3)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% averaged - left ROI (anterior) - t2star
subplot(8,2,2)
bar(subcode,corr(coords_averaged_left_ATL_ant_t2star_D3,U(:,3)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - t2star')
% averaged - left ROI (anterior) - tedana
subplot(8,2,4)
bar(subcode,corr(coords_averaged_left_ATL_ant_tedana_D3,U(:,3)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - tedana')
% averaged - left ROI (posterior) - t2star
subplot(8,2,6)
bar(subcode,corr(coords_averaged_left_ATL_pos_t2star_D3,U(:,3)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - t2star')
% averaged - left ROI (posterior) - tedana
subplot(8,2,8)
bar(subcode,corr(coords_averaged_left_ATL_pos_tedana_D3,U(:,3)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - tedana')
% averaged - right ROI (anterior) - t2star
subplot(8,2,10)
bar(subcode,corr(coords_averaged_right_ATL_ant_t2star_D3,U(:,3)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - t2star')
% averaged - right ROI (anterior) - tedana
subplot(8,2,12)
bar(subcode,corr(coords_averaged_right_ATL_ant_tedana_D3,U(:,3)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
hold on
% averaged - right ROI (posterior) - t2star
subplot(8,2,14)
bar(subcode,corr(coords_averaged_right_ATL_pos_t2star_D3,U(:,3)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Averaged - right ROI (posterior) - t2star')
% averaged - right ROI (posterior) - tedana
subplot(8,2,16)
bar(subcode,corr(coords_averaged_right_ATL_pos_tedana_D3,U(:,3)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
sgtitle('Correlations with predicted coordinates - all stimuli - dimension 3 - anterior and posterior halves of ROI')

% Dimension 1, animate only, anterior/posterior halves of ROI
figure;
% across run - left ROI (anterior) - t2star
subplot(8,2,1)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_t2star_D1(1:50,:),U(1:50,1)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Across run - left ROI (anterior) - t2star')
% across run - left ROI (anterior) - tedana
subplot(8,2,3)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_tedana_D1(1:50,:),U(1:50,1)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - left ROI (posterior) - t2star
subplot(8,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_t2star_D1(1:50,:),U(1:50,1)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Across run - left ROI (posterior) - t2star')
% across run - left ROI (posterior) - tedana
subplot(8,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_tedana_D1(1:50,:),U(1:50,1)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - right ROI (anterior) - t2star
subplot(8,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_t2star_D1(1:50,:),U(1:50,1)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ROI (anterior) - t2star')
% across run - right ROI (anterior) - tedana
subplot(8,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_tedana_D1(1:50,:),U(1:50,1)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% across run - right ROI (posterior) - t2star
subplot(8,2,13)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_t2star_D1(1:50,:),U(1:50,1)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Across run - right ROI (posterior) - t2star')
% across run - right ROI (posterior) - tedana
subplot(8,2,15)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_tedana_D1(1:50,:),U(1:50,1)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% averaged - left ROI (anterior) - t2star
subplot(8,2,2)
bar(subcode,corr(coords_averaged_left_ATL_ant_t2star_D1(1:50,:),U(1:50,1)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - t2star')
% averaged - left ROI (anterior) - tedana
subplot(8,2,4)
bar(subcode,corr(coords_averaged_left_ATL_ant_tedana_D1(1:50,:),U(1:50,1)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - tedana')
% averaged - left ROI (posterior) - t2star
subplot(8,2,6)
bar(subcode,corr(coords_averaged_left_ATL_pos_t2star_D1(1:50,:),U(1:50,1)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - t2star')
% averaged - left ROI (posterior) - tedana
subplot(8,2,8)
bar(subcode,corr(coords_averaged_left_ATL_pos_tedana_D1(1:50,:),U(1:50,1)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - tedana')
% averaged - right ROI (anterior) - t2star
subplot(8,2,10)
bar(subcode,corr(coords_averaged_right_ATL_ant_t2star_D1(1:50,:),U(1:50,1)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - t2star')
% averaged - right ROI (anterior) - tedana
subplot(8,2,12)
bar(subcode,corr(coords_averaged_right_ATL_ant_tedana_D1(1:50,:),U(1:50,1)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
hold on
% averaged - right ROI (posterior) - t2star
subplot(8,2,14)
bar(subcode,corr(coords_averaged_right_ATL_pos_t2star_D1(1:50,:),U(1:50,1)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Averaged - right ROI (posterior) - t2star')
% averaged - right ROI (posterior) - tedana
subplot(8,2,16)
bar(subcode,corr(coords_averaged_right_ATL_pos_tedana_D1(1:50,:),U(1:50,1)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
sgtitle('Correlations with predicted coordinates - animate only - dimension 1 - anterior and posterior halves of ROI')

% Dimension 2, animate only, anterior/posterior halves of ROI
figure;
% across run - left ROI (anterior) - t2star
subplot(8,2,1)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_t2star_D2(1:50,:),U(1:50,2)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Across run - left ROI (anterior) - t2star')
% across run - left ROI (anterior) - tedana
subplot(8,2,3)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_tedana_D2(1:50,:),U(1:50,2)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - left ROI (posterior) - t2star
subplot(8,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_t2star_D2(1:50,:),U(1:50,2)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Across run - left ROI (posterior) - t2star')
% across run - left ROI (posterior) - tedana
subplot(8,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_tedana_D2(1:50,:),U(1:50,2)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - right ROI (anterior) - t2star
subplot(8,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_t2star_D2(1:50,:),U(1:50,2)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ROI (anterior) - t2star')
% across run - right ROI (anterior) - tedana
subplot(8,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_tedana_D2(1:50,:),U(1:50,2)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% across run - right ROI (posterior) - t2star
subplot(8,2,13)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_t2star_D2(1:50,:),U(1:50,2)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Across run - right ROI (posterior) - t2star')
% across run - right ROI (posterior) - tedana
subplot(8,2,15)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_tedana_D2(1:50,:),U(1:50,2)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% averaged - left ROI (anterior) - t2star
subplot(8,2,2)
bar(subcode,corr(coords_averaged_left_ATL_ant_t2star_D2(1:50,:),U(1:50,2)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - t2star')
% averaged - left ROI (anterior) - tedana
subplot(8,2,4)
bar(subcode,corr(coords_averaged_left_ATL_ant_tedana_D2(1:50,:),U(1:50,2)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - tedana')
% averaged - left ROI (posterior) - t2star
subplot(8,2,6)
bar(subcode,corr(coords_averaged_left_ATL_pos_t2star_D2(1:50,:),U(1:50,2)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - t2star')
% averaged - left ROI (posterior) - tedana
subplot(8,2,8)
bar(subcode,corr(coords_averaged_left_ATL_pos_tedana_D2(1:50,:),U(1:50,2)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - tedana')
% averaged - right ROI (anterior) - t2star
subplot(8,2,10)
bar(subcode,corr(coords_averaged_right_ATL_ant_t2star_D2(1:50,:),U(1:50,2)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - t2star')
% averaged - right ROI (anterior) - tedana
subplot(8,2,12)
bar(subcode,corr(coords_averaged_right_ATL_ant_tedana_D2(1:50,:),U(1:50,2)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
hold on
% averaged - right ROI (posterior) - t2star
subplot(8,2,14)
bar(subcode,corr(coords_averaged_right_ATL_pos_t2star_D2(1:50,:),U(1:50,2)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Averaged - right ROI (posterior) - t2star')
% averaged - right ROI (posterior) - tedana
subplot(8,2,16)
bar(subcode,corr(coords_averaged_right_ATL_pos_tedana_D2(1:50,:),U(1:50,2)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
sgtitle('Correlations with predicted coordinates - animate only - dimension 2 - anterior and posterior halves of ROI')

% Dimension 3, animate only, anterior/posterior halves of ROI
figure;
% across run - left ROI (anterior) - t2star
subplot(8,2,1)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_t2star_D3(1:50,:),U(1:50,3)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Across run - left ROI (anterior) - t2star')
% across run - left ROI (anterior) - tedana
subplot(8,2,3)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_tedana_D3(1:50,:),U(1:50,3)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - left ROI (posterior) - t2star
subplot(8,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_t2star_D3(1:50,:),U(1:50,3)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Across run - left ROI (posterior) - t2star')
% across run - left ROI (posterior) - tedana
subplot(8,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_tedana_D3(1:50,:),U(1:50,3)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - right ROI (anterior) - t2star
subplot(8,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_t2star_D3(1:50,:),U(1:50,3)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ROI (anterior) - t2star')
% across run - right ROI (anterior) - tedana
subplot(8,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_tedana_D3(1:50,:),U(1:50,3)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% across run - right ROI (posterior) - t2star
subplot(8,2,13)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_t2star_D3(1:50,:),U(1:50,3)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Across run - right ROI (posterior) - t2star')
% across run - right ROI (posterior) - tedana
subplot(8,2,15)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_tedana_D3(1:50,:),U(1:50,3)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% averaged - left ROI (anterior) - t2star
subplot(8,2,2)
bar(subcode,corr(coords_averaged_left_ATL_ant_t2star_D3(1:50,:),U(1:50,3)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - t2star')
% averaged - left ROI (anterior) - tedana
subplot(8,2,4)
bar(subcode,corr(coords_averaged_left_ATL_ant_tedana_D3(1:50,:),U(1:50,3)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - tedana')
% averaged - left ROI (posterior) - t2star
subplot(8,2,6)
bar(subcode,corr(coords_averaged_left_ATL_pos_t2star_D3(1:50,:),U(1:50,3)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - t2star')
% averaged - left ROI (posterior) - tedana
subplot(8,2,8)
bar(subcode,corr(coords_averaged_left_ATL_pos_tedana_D3(1:50,:),U(1:50,3)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - tedana')
% averaged - right ROI (anterior) - t2star
subplot(8,2,10)
bar(subcode,corr(coords_averaged_right_ATL_ant_t2star_D3(1:50,:),U(1:50,3)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - t2star')
% averaged - right ROI (anterior) - tedana
subplot(8,2,12)
bar(subcode,corr(coords_averaged_right_ATL_ant_tedana_D3(1:50,:),U(1:50,3)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
hold on
% averaged - right ROI (posterior) - t2star
subplot(8,2,14)
bar(subcode,corr(coords_averaged_right_ATL_pos_t2star_D3(1:50,:),U(1:50,3)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Averaged - right ROI (posterior) - t2star')
% averaged - right ROI (posterior) - tedana
subplot(8,2,16)
bar(subcode,corr(coords_averaged_right_ATL_pos_tedana_D3(1:50,:),U(1:50,3)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
sgtitle('Correlations with predicted coordinates - animate only - dimension 3 - anterior and posterior halves of ROI')

% Dimension 1, inanimate only, anterior/posterior halves of ROI
figure;
% across run - left ROI (anterior) - t2star
subplot(8,2,1)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_t2star_D1(51:100,:),U(51:100,1)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Across run - left ROI (anterior) - t2star')
% across run - left ROI (anterior) - tedana
subplot(8,2,3)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_tedana_D1(51:100,:),U(51:100,1)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - left ROI (posterior) - t2star
subplot(8,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_t2star_D1(51:100,:),U(51:100,1)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Across run - left ROI (posterior) - t2star')
% across run - left ROI (posterior) - tedana
subplot(8,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_tedana_D1(51:100,:),U(51:100,1)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - right ROI (anterior) - t2star
subplot(8,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_t2star_D1(51:100,:),U(51:100,1)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ROI (anterior) - t2star')
% across run - right ROI (anterior) - tedana
subplot(8,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_tedana_D1(51:100,:),U(51:100,1)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% across run - right ROI (posterior) - t2star
subplot(8,2,13)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_t2star_D1(51:100,:),U(51:100,1)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Across run - right ROI (posterior) - t2star')
% across run - right ROI (posterior) - tedana
subplot(8,2,15)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_tedana_D1(51:100,:),U(51:100,1)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% averaged - left ROI (anterior) - t2star
subplot(8,2,2)
bar(subcode,corr(coords_averaged_left_ATL_ant_t2star_D1(51:100,:),U(51:100,1)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - t2star')
% averaged - left ROI (anterior) - tedana
subplot(8,2,4)
bar(subcode,corr(coords_averaged_left_ATL_ant_tedana_D1(51:100,:),U(51:100,1)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - tedana')
% averaged - left ROI (posterior) - t2star
subplot(8,2,6)
bar(subcode,corr(coords_averaged_left_ATL_pos_t2star_D1(51:100,:),U(51:100,1)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - t2star')
% averaged - left ROI (posterior) - tedana
subplot(8,2,8)
bar(subcode,corr(coords_averaged_left_ATL_pos_tedana_D1(51:100,:),U(51:100,1)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - tedana')
% averaged - right ROI (anterior) - t2star
subplot(8,2,10)
bar(subcode,corr(coords_averaged_right_ATL_ant_t2star_D1(51:100,:),U(51:100,1)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - t2star')
% averaged - right ROI (anterior) - tedana
subplot(8,2,12)
bar(subcode,corr(coords_averaged_right_ATL_ant_tedana_D1(51:100,:),U(51:100,1)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
hold on
% averaged - right ROI (posterior) - t2star
subplot(8,2,14)
bar(subcode,corr(coords_averaged_right_ATL_pos_t2star_D1(51:100,:),U(51:100,1)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Averaged - right ROI (posterior) - t2star')
% averaged - right ROI (posterior) - tedana
subplot(8,2,16)
bar(subcode,corr(coords_averaged_right_ATL_pos_tedana_D1(51:100,:),U(51:100,1)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
sgtitle('Correlations with predicted coordinates - inanimate only - dimension 1 - anterior and posterior halves of ROI')

% Dimension 2, inanimate only, anterior/posterior halves of ROI
figure;
% across run - left ROI (anterior) - t2star
subplot(8,2,1)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_t2star_D2(51:100,:),U(51:100,2)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Across run - left ROI (anterior) - t2star')
% across run - left ROI (anterior) - tedana
subplot(8,2,3)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_tedana_D2(51:100,:),U(51:100,2)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - left ROI (posterior) - t2star
subplot(8,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_t2star_D2(51:100,:),U(51:100,2)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Across run - left ROI (posterior) - t2star')
% across run - left ROI (posterior) - tedana
subplot(8,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_tedana_D2(51:100,:),U(51:100,2)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - right ROI (anterior) - t2star
subplot(8,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_t2star_D2(51:100,:),U(51:100,2)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ROI (anterior) - t2star')
% across run - right ROI (anterior) - tedana
subplot(8,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_tedana_D2(51:100,:),U(51:100,2)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% across run - right ROI (posterior) - t2star
subplot(8,2,13)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_t2star_D2(51:100,:),U(51:100,2)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Across run - right ROI (posterior) - t2star')
% across run - right ROI (posterior) - tedana
subplot(8,2,15)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_tedana_D2(51:100,:),U(51:100,2)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% averaged - left ROI (anterior) - t2star
subplot(8,2,2)
bar(subcode,corr(coords_averaged_left_ATL_ant_t2star_D2(51:100,:),U(51:100,2)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - t2star')
% averaged - left ROI (anterior) - tedana
subplot(8,2,4)
bar(subcode,corr(coords_averaged_left_ATL_ant_tedana_D2(51:100,:),U(51:100,2)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - tedana')
% averaged - left ROI (posterior) - t2star
subplot(8,2,6)
bar(subcode,corr(coords_averaged_left_ATL_pos_t2star_D2(51:100,:),U(51:100,2)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - t2star')
% averaged - left ROI (posterior) - tedana
subplot(8,2,8)
bar(subcode,corr(coords_averaged_left_ATL_pos_tedana_D2(51:100,:),U(51:100,2)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - tedana')
% averaged - right ROI (anterior) - t2star
subplot(8,2,10)
bar(subcode,corr(coords_averaged_right_ATL_ant_t2star_D2(51:100,:),U(51:100,2)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - t2star')
% averaged - right ROI (anterior) - tedana
subplot(8,2,12)
bar(subcode,corr(coords_averaged_right_ATL_ant_tedana_D2(51:100,:),U(51:100,2)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
hold on
% averaged - right ROI (posterior) - t2star
subplot(8,2,14)
bar(subcode,corr(coords_averaged_right_ATL_pos_t2star_D2(51:100,:),U(51:100,2)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Averaged - right ROI (posterior) - t2star')
% averaged - right ROI (posterior) - tedana
subplot(8,2,16)
bar(subcode,corr(coords_averaged_right_ATL_pos_tedana_D2(51:100,:),U(51:100,2)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
sgtitle('Correlations with predicted coordinates - inanimate only - dimension 2 - anterior and posterior halves of ROI')

% Dimension 3, inanimate only, anterior/posterior halves of ROI
figure;
% across run - left ROI (anterior) - t2star
subplot(8,2,1)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_t2star_D3(51:100,:),U(51:100,3)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Across run - left ROI (anterior) - t2star')
% across run - left ROI (anterior) - tedana
subplot(8,2,3)
bar(subcode,corr(coords_acrossRun_left_ATL_ant_tedana_D3(51:100,:),U(51:100,3)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - left ROI (posterior) - t2star
subplot(8,2,5)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_t2star_D3(51:100,:),U(51:100,3)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Across run - left ROI (posterior) - t2star')
% across run - left ROI (posterior) - tedana
subplot(8,2,7)
bar(subcode,corr(coords_acrossRun_left_ATL_pos_tedana_D3(51:100,:),U(51:100,3)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - left ROI (anterior) - tedana')
% across run - right ROI (anterior) - t2star
subplot(8,2,9)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_t2star_D3(51:100,:),U(51:100,3)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Across run - right ROI (anterior) - t2star')
% across run - right ROI (anterior) - tedana
subplot(8,2,11)
bar(subcode,corr(coords_acrossRun_right_ATL_ant_tedana_D3(51:100,:),U(51:100,3)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% across run - right ROI (posterior) - t2star
subplot(8,2,13)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_t2star_D3(51:100,:),U(51:100,3)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Across run - right ROI (posterior) - t2star')
% across run - right ROI (posterior) - tedana
subplot(8,2,15)
bar(subcode,corr(coords_acrossRun_right_ATL_pos_tedana_D3(51:100,:),U(51:100,3)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Across run - right ROI (anterior) - tedana')
% averaged - left ROI (anterior) - t2star
subplot(8,2,2)
bar(subcode,corr(coords_averaged_left_ATL_ant_t2star_D3(51:100,:),U(51:100,3)),'FaceColor','#D95319')
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - t2star')
% averaged - left ROI (anterior) - tedana
subplot(8,2,4)
bar(subcode,corr(coords_averaged_left_ATL_ant_tedana_D3(51:100,:),U(51:100,3)),'FaceColor','#D95319','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (anterior) - tedana')
% averaged - left ROI (posterior) - t2star
subplot(8,2,6)
bar(subcode,corr(coords_averaged_left_ATL_pos_t2star_D3(51:100,:),U(51:100,3)),'FaceColor','#EDB120')
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - t2star')
% averaged - left ROI (posterior) - tedana
subplot(8,2,8)
bar(subcode,corr(coords_averaged_left_ATL_pos_tedana_D3(51:100,:),U(51:100,3)),'FaceColor','#EDB120','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - left ROI (posterior) - tedana')
% averaged - right ROI (anterior) - t2star
subplot(8,2,10)
bar(subcode,corr(coords_averaged_right_ATL_ant_t2star_D3(51:100,:),U(51:100,3)),'FaceColor','#0072BD')
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - t2star')
% averaged - right ROI (anterior) - tedana
subplot(8,2,12)
bar(subcode,corr(coords_averaged_right_ATL_ant_tedana_D3(51:100,:),U(51:100,3)),'FaceColor','#0072BD','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
hold on
% averaged - right ROI (posterior) - t2star
subplot(8,2,14)
bar(subcode,corr(coords_averaged_right_ATL_pos_t2star_D3(51:100,:),U(51:100,3)),'FaceColor','#4DBEEE')
ylim([-0.1,1])
title('Averaged - right ROI (posterior) - t2star')
% averaged - right ROI (posterior) - tedana
subplot(8,2,16)
bar(subcode,corr(coords_averaged_right_ATL_pos_tedana_D3(51:100,:),U(51:100,3)),'FaceColor','#4DBEEE','FaceAlpha',0.25)
ylim([-0.1,1])
title('Averaged - right ROI (anterior) - tedana')
sgtitle('Correlations with predicted coordinates - inanimate only - dimension 3 - anterior and posterior halves of ROI')

% group results, all stimuli, subdivided ROIs
groupsubcorrall = zeros(size(subcode,2),48);
groupsubcorrall(:,1) = corr(coords_acrossRun_left_ATL_ant_t2star_D1,U(:,1));
groupsubcorrall(:,2) = corr(coords_acrossRun_left_ATL_ant_t2star_D2,U(:,2));
groupsubcorrall(:,3) = corr(coords_acrossRun_left_ATL_ant_t2star_D3,U(:,3));
groupsubcorrall(:,4) = corr(coords_acrossRun_left_ATL_ant_tedana_D1,U(:,1));
groupsubcorrall(:,5) = corr(coords_acrossRun_left_ATL_ant_tedana_D2,U(:,2));
groupsubcorrall(:,6) = corr(coords_acrossRun_left_ATL_ant_tedana_D3,U(:,3));
groupsubcorrall(:,7) = corr(coords_acrossRun_left_ATL_pos_t2star_D1,U(:,1));
groupsubcorrall(:,8) = corr(coords_acrossRun_left_ATL_pos_t2star_D2,U(:,2));
groupsubcorrall(:,9) = corr(coords_acrossRun_left_ATL_pos_t2star_D3,U(:,3));
groupsubcorrall(:,10) = corr(coords_acrossRun_left_ATL_pos_tedana_D1,U(:,1));
groupsubcorrall(:,11) = corr(coords_acrossRun_left_ATL_pos_tedana_D2,U(:,2));
groupsubcorrall(:,12) = corr(coords_acrossRun_left_ATL_pos_tedana_D3,U(:,3));
groupsubcorrall(:,13) = corr(coords_acrossRun_right_ATL_ant_t2star_D1,U(:,1));
groupsubcorrall(:,14) = corr(coords_acrossRun_right_ATL_ant_t2star_D2,U(:,2));
groupsubcorrall(:,15) = corr(coords_acrossRun_right_ATL_ant_t2star_D3,U(:,3));
groupsubcorrall(:,16) = corr(coords_acrossRun_right_ATL_ant_tedana_D1,U(:,1));
groupsubcorrall(:,17) = corr(coords_acrossRun_right_ATL_ant_tedana_D2,U(:,2));
groupsubcorrall(:,18) = corr(coords_acrossRun_right_ATL_ant_tedana_D3,U(:,3));
groupsubcorrall(:,19) = corr(coords_acrossRun_right_ATL_pos_t2star_D1,U(:,1));
groupsubcorrall(:,20) = corr(coords_acrossRun_right_ATL_pos_t2star_D2,U(:,2));
groupsubcorrall(:,21) = corr(coords_acrossRun_right_ATL_pos_t2star_D3,U(:,3));
groupsubcorrall(:,22) = corr(coords_acrossRun_right_ATL_pos_tedana_D1,U(:,1));
groupsubcorrall(:,23) = corr(coords_acrossRun_right_ATL_pos_tedana_D2,U(:,2));
groupsubcorrall(:,24) = corr(coords_acrossRun_right_ATL_pos_tedana_D3,U(:,3));
groupsubcorrall(:,25) = corr(coords_averaged_left_ATL_ant_t2star_D1,U(:,1));
groupsubcorrall(:,26) = corr(coords_averaged_left_ATL_ant_t2star_D2,U(:,2));
groupsubcorrall(:,27) = corr(coords_averaged_left_ATL_ant_t2star_D3,U(:,3));
groupsubcorrall(:,28) = corr(coords_averaged_left_ATL_ant_tedana_D1,U(:,1));
groupsubcorrall(:,29) = corr(coords_averaged_left_ATL_ant_tedana_D2,U(:,2));
groupsubcorrall(:,30) = corr(coords_averaged_left_ATL_ant_tedana_D3,U(:,3));
groupsubcorrall(:,31) = corr(coords_averaged_left_ATL_pos_t2star_D1,U(:,1));
groupsubcorrall(:,32) = corr(coords_averaged_left_ATL_pos_t2star_D2,U(:,2));
groupsubcorrall(:,33) = corr(coords_averaged_left_ATL_pos_t2star_D3,U(:,3));
groupsubcorrall(:,34) = corr(coords_averaged_left_ATL_pos_tedana_D1,U(:,1));
groupsubcorrall(:,35) = corr(coords_averaged_left_ATL_pos_tedana_D2,U(:,2));
groupsubcorrall(:,36) = corr(coords_averaged_left_ATL_pos_tedana_D3,U(:,3));
groupsubcorrall(:,37) = corr(coords_averaged_right_ATL_ant_t2star_D1,U(:,1));
groupsubcorrall(:,38) = corr(coords_averaged_right_ATL_ant_t2star_D2,U(:,2));
groupsubcorrall(:,39) = corr(coords_averaged_right_ATL_ant_t2star_D3,U(:,3));
groupsubcorrall(:,40) = corr(coords_averaged_right_ATL_ant_tedana_D1,U(:,1));
groupsubcorrall(:,41) = corr(coords_averaged_right_ATL_ant_tedana_D2,U(:,2));
groupsubcorrall(:,42) = corr(coords_averaged_right_ATL_ant_tedana_D3,U(:,3));
groupsubcorrall(:,43) = corr(coords_averaged_right_ATL_pos_t2star_D1,U(:,1));
groupsubcorrall(:,44) = corr(coords_averaged_right_ATL_pos_t2star_D2,U(:,2));
groupsubcorrall(:,45) = corr(coords_averaged_right_ATL_pos_t2star_D3,U(:,3));
groupsubcorrall(:,46) = corr(coords_averaged_right_ATL_pos_tedana_D1,U(:,1));
groupsubcorrall(:,47) = corr(coords_averaged_right_ATL_pos_tedana_D2,U(:,2));
groupsubcorrall(:,48) = corr(coords_averaged_right_ATL_pos_tedana_D3,U(:,3));
colnames = {'AR-lATL-ant-t2star-D1','AR-lATL-ant-t2star-D2','AR-lATL-ant-t2star-D3','AR-lATL-ant-tedana-D1','AR-lATL-ant-tedana-D2','AR-lATL-ant-tedana-D3','AR-lATL-pos-t2star-D1','AR-lATL-pos-t2star-D2','AR-lATL-pos-t2star-D3','AR-lATL-pos-tedana-D1','AR-lATL-pos-tedana-D2','AR-lATL-pos-tedana-D3','AR-rATL-ant-t2star-D1','AR-rATL-ant-t2star-D2','AR-rATL-ant-t2star-D3','AR-rATL-ant-tedana-D1','AR-rATL-ant-tedana-D2','AR-rATL-ant-tedana-D3','AR-rATL-pos-t2star-D1','AR-rATL-pos-t2star-D2','AR-rATL-pos-t2star-D3','AR-rATL-pos-tedana-D1','AR-rATL-pos-tedana-D2','AR-rATL-pos-tedana-D3','av-lATL-ant-t2star-D1','av-lATL-ant-t2star-D2','av-lATL-ant-t2star-D3','av-lATL-ant-tedana-D1','av-lATL-ant-tedana-D2','av-lATL-ant-tedana-D3','av-lATL-pos-t2star-D1','av-lATL-pos-t2star-D2','av-lATL-pos-t2star-D3','av-lATL-pos-tedana-D1','av-lATL-pos-tedana-D2','av-lATL-pos-tedana-D3','av-rATL-ant-t2star-D1','av-rATL-ant-t2star-D2','av-rATL-ant-t2star-D3','av-rATL-ant-tedana-D1','av-rATL-ant-tedana-D2','av-rATL-ant-tedana-D3','av-rATL-pos-t2star-D1','av-rATL-pos-t2star-D2','av-rATL-pos-t2star-D3','av-rATL-pos-tedana-D1','av-rATL-pos-tedana-D2','av-rATL-pos-tedana-D3',};

figure;
b = bar(colnames,mean(groupsubcorrall),'FaceColor','flat');
% normalised RGB triplets
colours = [0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330];
b.CData = colours;
ylim([-0.1,1])
title('Group results - correlations with predicted coordinates')
hold on
confint = 1.96*(std(groupsubcorrall)/sqrt(size(groupsubcorrall,1)));
errorbar(1:size(groupsubcorrall,2),mean(groupsubcorrall,1),confint,'LineStyle','none','Color','black');

% group results, animate only, subdivided ROIs
groupsubcorranimate = zeros(size(subcode,2),48);
groupsubcorranimate(:,1) = corr(coords_acrossRun_left_ATL_ant_t2star_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,2) = corr(coords_acrossRun_left_ATL_ant_t2star_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,3) = corr(coords_acrossRun_left_ATL_ant_t2star_D3(1:50,:),U(1:50,3));
groupsubcorranimate(:,4) = corr(coords_acrossRun_left_ATL_ant_tedana_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,5) = corr(coords_acrossRun_left_ATL_ant_tedana_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,6) = corr(coords_acrossRun_left_ATL_ant_tedana_D3(1:50,:),U(1:50,3));
groupsubcorranimate(:,7) = corr(coords_acrossRun_left_ATL_pos_t2star_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,8) = corr(coords_acrossRun_left_ATL_pos_t2star_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,9) = corr(coords_acrossRun_left_ATL_pos_t2star_D3(1:50,:),U(1:50,3));
groupsubcorranimate(:,10) = corr(coords_acrossRun_left_ATL_pos_tedana_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,11) = corr(coords_acrossRun_left_ATL_pos_tedana_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,12) = corr(coords_acrossRun_left_ATL_pos_tedana_D3(1:50,:),U(1:50,3));
groupsubcorranimate(:,13) = corr(coords_acrossRun_right_ATL_ant_t2star_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,14) = corr(coords_acrossRun_right_ATL_ant_t2star_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,15) = corr(coords_acrossRun_right_ATL_ant_t2star_D3(1:50,:),U(1:50,3));
groupsubcorranimate(:,16) = corr(coords_acrossRun_right_ATL_ant_tedana_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,17) = corr(coords_acrossRun_right_ATL_ant_tedana_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,18) = corr(coords_acrossRun_right_ATL_ant_tedana_D3(1:50,:),U(1:50,3));
groupsubcorranimate(:,19) = corr(coords_acrossRun_right_ATL_pos_t2star_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,20) = corr(coords_acrossRun_right_ATL_pos_t2star_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,21) = corr(coords_acrossRun_right_ATL_pos_t2star_D3(1:50,:),U(1:50,3));
groupsubcorranimate(:,22) = corr(coords_acrossRun_right_ATL_pos_tedana_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,23) = corr(coords_acrossRun_right_ATL_pos_tedana_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,24) = corr(coords_acrossRun_right_ATL_pos_tedana_D3(1:50,:),U(1:50,3));
groupsubcorranimate(:,25) = corr(coords_averaged_left_ATL_ant_t2star_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,26) = corr(coords_averaged_left_ATL_ant_t2star_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,27) = corr(coords_averaged_left_ATL_ant_t2star_D3(1:50,:),U(1:50,3));
groupsubcorranimate(:,28) = corr(coords_averaged_left_ATL_ant_tedana_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,29) = corr(coords_averaged_left_ATL_ant_tedana_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,30) = corr(coords_averaged_left_ATL_ant_tedana_D3(1:50,:),U(1:50,3));
groupsubcorranimate(:,31) = corr(coords_averaged_left_ATL_pos_t2star_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,32) = corr(coords_averaged_left_ATL_pos_t2star_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,33) = corr(coords_averaged_left_ATL_pos_t2star_D3(1:50,:),U(1:50,3));
groupsubcorranimate(:,34) = corr(coords_averaged_left_ATL_pos_tedana_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,35) = corr(coords_averaged_left_ATL_pos_tedana_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,36) = corr(coords_averaged_left_ATL_pos_tedana_D3(1:50,:),U(1:50,3));
groupsubcorranimate(:,37) = corr(coords_averaged_right_ATL_ant_t2star_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,38) = corr(coords_averaged_right_ATL_ant_t2star_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,39) = corr(coords_averaged_right_ATL_ant_t2star_D3(1:50,:),U(1:50,3));
groupsubcorranimate(:,40) = corr(coords_averaged_right_ATL_ant_tedana_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,41) = corr(coords_averaged_right_ATL_ant_tedana_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,42) = corr(coords_averaged_right_ATL_ant_tedana_D3(1:50,:),U(1:50,3));
groupsubcorranimate(:,43) = corr(coords_averaged_right_ATL_pos_t2star_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,44) = corr(coords_averaged_right_ATL_pos_t2star_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,45) = corr(coords_averaged_right_ATL_pos_t2star_D3(1:50,:),U(1:50,3));
groupsubcorranimate(:,46) = corr(coords_averaged_right_ATL_pos_tedana_D1(1:50,:),U(1:50,1));
groupsubcorranimate(:,47) = corr(coords_averaged_right_ATL_pos_tedana_D2(1:50,:),U(1:50,2));
groupsubcorranimate(:,48) = corr(coords_averaged_right_ATL_pos_tedana_D3(1:50,:),U(1:50,3));
colnames = {'AR-lATL-ant-t2star-D1','AR-lATL-ant-t2star-D2','AR-lATL-ant-t2star-D3','AR-lATL-ant-tedana-D1','AR-lATL-ant-tedana-D2','AR-lATL-ant-tedana-D3','AR-lATL-pos-t2star-D1','AR-lATL-pos-t2star-D2','AR-lATL-pos-t2star-D3','AR-lATL-pos-tedana-D1','AR-lATL-pos-tedana-D2','AR-lATL-pos-tedana-D3','AR-rATL-ant-t2star-D1','AR-rATL-ant-t2star-D2','AR-rATL-ant-t2star-D3','AR-rATL-ant-tedana-D1','AR-rATL-ant-tedana-D2','AR-rATL-ant-tedana-D3','AR-rATL-pos-t2star-D1','AR-rATL-pos-t2star-D2','AR-rATL-pos-t2star-D3','AR-rATL-pos-tedana-D1','AR-rATL-pos-tedana-D2','AR-rATL-pos-tedana-D3','av-lATL-ant-t2star-D1','av-lATL-ant-t2star-D2','av-lATL-ant-t2star-D3','av-lATL-ant-tedana-D1','av-lATL-ant-tedana-D2','av-lATL-ant-tedana-D3','av-lATL-pos-t2star-D1','av-lATL-pos-t2star-D2','av-lATL-pos-t2star-D3','av-lATL-pos-tedana-D1','av-lATL-pos-tedana-D2','av-lATL-pos-tedana-D3','av-rATL-ant-t2star-D1','av-rATL-ant-t2star-D2','av-rATL-ant-t2star-D3','av-rATL-ant-tedana-D1','av-rATL-ant-tedana-D2','av-rATL-ant-tedana-D3','av-rATL-pos-t2star-D1','av-rATL-pos-t2star-D2','av-rATL-pos-t2star-D3','av-rATL-pos-tedana-D1','av-rATL-pos-tedana-D2','av-rATL-pos-tedana-D3',};

figure;
b = bar(colnames,mean(groupsubcorranimate),'FaceColor','flat');
% normalised RGB triplets
colours = [0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330];
b.CData = colours;
ylim([-0.1,1])
title('Group results - animate only - correlations with predicted coordinates')
hold on
confint = 1.96*(std(groupsubcorranimate)/sqrt(size(groupsubcorranimate,1)));
errorbar(1:size(groupsubcorranimate,2),mean(groupsubcorranimate,1),confint,'LineStyle','none','Color','black');

% group results, inanimate only, subdivided ROIs
groupsubcorrinanimate = zeros(size(subcode,2),48);
groupsubcorrinanimate(:,1) = corr(coords_acrossRun_left_ATL_ant_t2star_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,2) = corr(coords_acrossRun_left_ATL_ant_t2star_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,3) = corr(coords_acrossRun_left_ATL_ant_t2star_D3(51:100,:),U(51:100,3));
groupsubcorrinanimate(:,4) = corr(coords_acrossRun_left_ATL_ant_tedana_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,5) = corr(coords_acrossRun_left_ATL_ant_tedana_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,6) = corr(coords_acrossRun_left_ATL_ant_tedana_D3(51:100,:),U(51:100,3));
groupsubcorrinanimate(:,7) = corr(coords_acrossRun_left_ATL_pos_t2star_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,8) = corr(coords_acrossRun_left_ATL_pos_t2star_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,9) = corr(coords_acrossRun_left_ATL_pos_t2star_D3(51:100,:),U(51:100,3));
groupsubcorrinanimate(:,10) = corr(coords_acrossRun_left_ATL_pos_tedana_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,11) = corr(coords_acrossRun_left_ATL_pos_tedana_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,12) = corr(coords_acrossRun_left_ATL_pos_tedana_D3(51:100,:),U(51:100,3));
groupsubcorrinanimate(:,13) = corr(coords_acrossRun_right_ATL_ant_t2star_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,14) = corr(coords_acrossRun_right_ATL_ant_t2star_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,15) = corr(coords_acrossRun_right_ATL_ant_t2star_D3(51:100,:),U(51:100,3));
groupsubcorrinanimate(:,16) = corr(coords_acrossRun_right_ATL_ant_tedana_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,17) = corr(coords_acrossRun_right_ATL_ant_tedana_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,18) = corr(coords_acrossRun_right_ATL_ant_tedana_D3(51:100,:),U(51:100,3));
groupsubcorrinanimate(:,19) = corr(coords_acrossRun_right_ATL_pos_t2star_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,20) = corr(coords_acrossRun_right_ATL_pos_t2star_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,21) = corr(coords_acrossRun_right_ATL_pos_t2star_D3(51:100,:),U(51:100,3));
groupsubcorrinanimate(:,22) = corr(coords_acrossRun_right_ATL_pos_tedana_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,23) = corr(coords_acrossRun_right_ATL_pos_tedana_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,24) = corr(coords_acrossRun_right_ATL_pos_tedana_D3(51:100,:),U(51:100,3));
groupsubcorrinanimate(:,25) = corr(coords_averaged_left_ATL_ant_t2star_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,26) = corr(coords_averaged_left_ATL_ant_t2star_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,27) = corr(coords_averaged_left_ATL_ant_t2star_D3(51:100,:),U(51:100,3));
groupsubcorrinanimate(:,28) = corr(coords_averaged_left_ATL_ant_tedana_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,29) = corr(coords_averaged_left_ATL_ant_tedana_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,30) = corr(coords_averaged_left_ATL_ant_tedana_D3(51:100,:),U(51:100,3));
groupsubcorrinanimate(:,31) = corr(coords_averaged_left_ATL_pos_t2star_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,32) = corr(coords_averaged_left_ATL_pos_t2star_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,33) = corr(coords_averaged_left_ATL_pos_t2star_D3(51:100,:),U(51:100,3));
groupsubcorrinanimate(:,34) = corr(coords_averaged_left_ATL_pos_tedana_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,35) = corr(coords_averaged_left_ATL_pos_tedana_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,36) = corr(coords_averaged_left_ATL_pos_tedana_D3(51:100,:),U(51:100,3));
groupsubcorrinanimate(:,37) = corr(coords_averaged_right_ATL_ant_t2star_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,38) = corr(coords_averaged_right_ATL_ant_t2star_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,39) = corr(coords_averaged_right_ATL_ant_t2star_D3(51:100,:),U(51:100,3));
groupsubcorrinanimate(:,40) = corr(coords_averaged_right_ATL_ant_tedana_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,41) = corr(coords_averaged_right_ATL_ant_tedana_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,42) = corr(coords_averaged_right_ATL_ant_tedana_D3(51:100,:),U(51:100,3));
groupsubcorrinanimate(:,43) = corr(coords_averaged_right_ATL_pos_t2star_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,44) = corr(coords_averaged_right_ATL_pos_t2star_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,45) = corr(coords_averaged_right_ATL_pos_t2star_D3(51:100,:),U(51:100,3));
groupsubcorrinanimate(:,46) = corr(coords_averaged_right_ATL_pos_tedana_D1(51:100,:),U(51:100,1));
groupsubcorrinanimate(:,47) = corr(coords_averaged_right_ATL_pos_tedana_D2(51:100,:),U(51:100,2));
groupsubcorrinanimate(:,48) = corr(coords_averaged_right_ATL_pos_tedana_D3(51:100,:),U(51:100,3));
colnames = {'AR-lATL-ant-t2star-D1','AR-lATL-ant-t2star-D2','AR-lATL-ant-t2star-D3','AR-lATL-ant-tedana-D1','AR-lATL-ant-tedana-D2','AR-lATL-ant-tedana-D3','AR-lATL-pos-t2star-D1','AR-lATL-pos-t2star-D2','AR-lATL-pos-t2star-D3','AR-lATL-pos-tedana-D1','AR-lATL-pos-tedana-D2','AR-lATL-pos-tedana-D3','AR-rATL-ant-t2star-D1','AR-rATL-ant-t2star-D2','AR-rATL-ant-t2star-D3','AR-rATL-ant-tedana-D1','AR-rATL-ant-tedana-D2','AR-rATL-ant-tedana-D3','AR-rATL-pos-t2star-D1','AR-rATL-pos-t2star-D2','AR-rATL-pos-t2star-D3','AR-rATL-pos-tedana-D1','AR-rATL-pos-tedana-D2','AR-rATL-pos-tedana-D3','av-lATL-ant-t2star-D1','av-lATL-ant-t2star-D2','av-lATL-ant-t2star-D3','av-lATL-ant-tedana-D1','av-lATL-ant-tedana-D2','av-lATL-ant-tedana-D3','av-lATL-pos-t2star-D1','av-lATL-pos-t2star-D2','av-lATL-pos-t2star-D3','av-lATL-pos-tedana-D1','av-lATL-pos-tedana-D2','av-lATL-pos-tedana-D3','av-rATL-ant-t2star-D1','av-rATL-ant-t2star-D2','av-rATL-ant-t2star-D3','av-rATL-ant-tedana-D1','av-rATL-ant-tedana-D2','av-rATL-ant-tedana-D3','av-rATL-pos-t2star-D1','av-rATL-pos-t2star-D2','av-rATL-pos-t2star-D3','av-rATL-pos-tedana-D1','av-rATL-pos-tedana-D2','av-rATL-pos-tedana-D3',};

figure;
b = bar(colnames,mean(groupsubcorrinanimate),'FaceColor','flat');
% normalised RGB triplets
colours = [0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0.9290,0.6940,0.1250;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0,0.4470,0.7410;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330;0.3010,0.7450,0.9330];
b.CData = colours;
ylim([-0.1,1])
title('Group results - inanimate only - correlations with predicted coordinates')
hold on
confint = 1.96*(std(groupsubcorrinanimate)/sqrt(size(groupsubcorrinanimate,1)));
errorbar(1:size(groupsubcorrinanimate,2),mean(groupsubcorrinanimate,1),confint,'LineStyle','none','Color','black');
