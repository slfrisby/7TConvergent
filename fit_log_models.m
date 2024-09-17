function fit_log_models(ids)
% by Saskia. Fits LASSO regression models by calling get_log_hoerr. 

% Current functionality:
% - fits models to both t2star and tedana data
% - fits models in two ways:
% - - average betas across 4 runs and do standard 10-fold cross-validation
% - - don't average across runs. Draw the test set from run 1 and the
%     training sets from runs 2, 3 and 4 (individual items are either in
%     the training set or the test set). Repeat with run 2 (same training
%     and test sets), then run 3, then run 4 to get 40-fold
%     cross-validation. 
% - fits models to whole brain and to 2 ECoG-inspired ATL ROIs (left and
%   right).


    root = ['/imaging/projects/cbu/wbic-p00591-DAISY/main/derivatives'];
    cd([root]);
    
    subcode{1} = ['sub-',ids] 
    % data = {'t2star','tedana'};
    data= {'tedana'};
    
    for s = 1:size(subcode,2)
    
        % for each data type
        for d = 1:length(data)

            % make output directory
            mkdir([root,'/LASSO/log/',data{d},'/',subcode{s}]);
    
            % load data matrix X
            disp(['Loading ',data{d},' data...'])
            load([root,'/cox/',subcode{s},'/',subcode{s},'_rec-',data{d},'_X.mat']);
            % load metadata, which contains ROIs set as filters
            load([root,'/cox/',data{d},'_full_metadata.mat']);

            % get position of participant's metadata within metadata
            % variable
            tmp = str2num(erase(subcode{1},'sub-'));
            subs = [];
            for i = 1:size(metadata,2)
                subs = [subs metadata(i).subject];
            end
            subidx = find(subs==tmp)
    
            % % 1. ACROSS-RUN
            % 
            % % set Y (vector of categories)
            % Y = repmat([zeros(50,1);ones(50,1)],4,1);
            % 
            % % 1A. Whole brain
            % disp('Decoding from the whole brain (across runs, 40-fold)...')
            % acrossRun_whole_brain = get_log_hoerr(X,Y,'cvind',metadata(subidx).cvind(:,1));
            % 
            % % 1B. Left ATL
            % leftATL = X(:,metadata(subidx).filters(5).filter);
            % disp('Decoding from the left ATL (across runs, 40-fold)...')
            % acrossRun_left_ATL = get_log_hoerr(leftATL,Y,'cvind',metadata(subidx).cvind(:,1));
            % 
            % % 1C. Right ATL
            % rightATL = X(:,metadata(subidx).filters(6).filter);
            % disp('Decoding from the right ATL (across runs, 40-fold)...')
            % acrossRun_right_ATL = get_log_hoerr(rightATL,Y,'cvind',metadata(subidx).cvind(:,1));
            % 
            % % 1D. Left ATL (anterior)
            % leftATLant = X(:,metadata(subidx).filters(7).filter);
            % disp('Decoding from the left ATL - anterior part only (across runs, 40-fold)...')
            % acrossRun_left_ATL_ant = get_log_hoerr(leftATLant,Y,'cvind',metadata(subidx).cvind(:,1));
            % 
            % % 1E. Left ATL (posterior)
            % leftATLpos = X(:,metadata(subidx).filters(8).filter);
            % disp('Decoding from the left ATL - posterior part only (across runs, 40-fold)...')
            % acrossRun_left_ATL_pos = get_log_hoerr(leftATLpos,Y,'cvind',metadata(subidx).cvind(:,1));
            % 
            % % 1F. Right ATL (anterior)
            % rightATLant = X(:,metadata(subidx).filters(9).filter);
            % disp('Decoding from the right ATL - anterior part only (across runs, 40-fold)...')
            % acrossRun_right_ATL_ant = get_log_hoerr(rightATLant,Y,'cvind',metadata(subidx).cvind(:,1));
            % 
            % % 1G. Right ATL (posterior)
            % rightATLpos = X(:,metadata(subidx).filters(10).filter);
            % disp('Decoding from the right ATL - posterior part only (across runs, 40-fold)...')
            % acrossRun_right_ATL_pos = get_log_hoerr(rightATLpos,Y,'cvind',metadata(subidx).cvind(:,1));

            % 2. AVERAGED

            % average X
            X = cat(3,X(1:100,:),X(101:200,:),X(201:300,:),X(301:400,:));
            X = mean(X,3);
            % load metadata, which contains ROIs set as filters
            load([root,'/cox/',data{d},'_averaged_metadata.mat']);
            % load PERMUTATION_STRUCT
            load([root,'/cox/',data{d},'_averaged_PERMUTATION_STRUCT.mat']);
    
            % set Y
            Y = [zeros(50,1);ones(50,1)];

            % 2A. Whole brain
            disp('Decoding from the whole brain (averaged across runs, 10-fold)...')
            averaged_whole_brain = get_log_hoerr(X,Y,'acrossRun',0,'cvind',metadata(subidx).cvind(:,1),'permstruct',PERMUTATION_STRUCT(subidx).permutation_index);

            % 2B. Left ATL
            leftATL = X(:,metadata(subidx).filters(5).filter);
            disp('Decoding from the left ATL (averaged across runs, 10-fold)...')
            averaged_left_ATL = get_log_hoerr(leftATL,Y,'acrossRun',0,'cvind',metadata(subidx).cvind(:,1),'permstruct',PERMUTATION_STRUCT(subidx).permutation_index);

            % 2C. Right ATL
            rightATL = X(:,metadata(subidx).filters(6).filter);
            disp('Decoding from the right ATL (averaged across runs, 10-fold)...')
            averaged_right_ATL = get_log_hoerr(rightATL,Y,'acrossRun',0,'cvind',metadata(subidx).cvind(:,1),'permstruct',PERMUTATION_STRUCT(subidx).permutation_index);

            % 2D. Left ATL (anterior)
            leftATLant = X(:,metadata(subidx).filters(7).filter);
            disp('Decoding from the left ATL - anterior part only (averaged across runs, 10-fold)...')
            averaged_left_ATL_ant = get_log_hoerr(leftATLant,Y,'acrossRun',0,'cvind',metadata(subidx).cvind(:,1),'permstruct',PERMUTATION_STRUCT(subidx).permutation_index);

            % 2E. Left ATL (posterior)
            leftATLpos = X(:,metadata(subidx).filters(8).filter);
            disp('Decoding from the left ATL - posterior part only (averaged across runs, 10-fold)...')
            averaged_left_ATL_pos = get_log_hoerr(leftATLpos,Y,'acrossRun',0,'cvind',metadata(subidx).cvind(:,1),'permstruct',PERMUTATION_STRUCT(subidx).permutation_index);

            % 2F. Right ATL (anterior)
            rightATLant = X(:,metadata(subidx).filters(9).filter);
            disp('Decoding from the right ATL - anterior part only (averaged across runs, 10-fold)...')
            averaged_right_ATL_ant = get_log_hoerr(rightATLant,Y,'acrossRun',0,'cvind',metadata(subidx).cvind(:,1),'permstruct',PERMUTATION_STRUCT(subidx).permutation_index);

            % 2G. Right ATL (posterior)
            rightATLpos = X(:,metadata(subidx).filters(10).filter);
            disp('Decoding from the right ATL - posterior part only (averaged across runs, 10-fold)...')
            averaged_right_ATL_pos = get_log_hoerr(rightATLpos,Y,'acrossRun',0,'cvind',metadata(subidx).cvind(:,1),'permstruct',PERMUTATION_STRUCT(subidx).permutation_index);
    
            % save everything in a single .mat file. Just save everything
            % (too many variables to faff about with individually!)
            save([root,'/LASSO/log/',data{d},'/',subcode{s},'/results.mat']);
        end
    end
end







