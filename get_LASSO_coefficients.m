function get_LASSO_coefficients(ids)
    % by Saskia. Gets coefficients for whole-brain LASSO models, linear and
    % logistic. 

    addpath(genpath('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/WISC_MVPA/'));

    root = ['/imaging/projects/cbu/wbic-p00591-DAISY/main/'];
    cd([root]);
    
    subcode{1} = ['sub-',ids] 
    % data = {'t2star','tedana'};
    data= {'tedana'};

    % decompose similarity matrix into 3 singular values
    load('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/dilkina_norms.mat');
    [C,z] = embed_similarity_matrix(dilkina_norms,3);
    U = rescale_embedding(C,z);

    for s = 1:size(subcode,2)
    
        % for each data type
        for d = 1:length(data)

            % make output directory
            mkdir([root,'/work/',subcode{s},'/coefficients/mat/',data{d},'/']);

            % load data matrix X
            disp(['Loading ',data{d},' data...'])
            load([root,'/derivatives/cox/',subcode{s},'/',subcode{s},'_rec-',data{d},'_X.mat']);
            % average X
            X = cat(3,X(1:100,:),X(101:200,:),X(201:300,:),X(301:400,:));
            X = mean(X,3);
            % load metadata, which contains ROIs set as filters
            load([root,'/derivatives/cox/',data{d},'_averaged_metadata.mat']);
            % load PERMUTATION_STRUCT
            load([root,'/derivatives/cox/',data{d},'_averaged_PERMUTATION_STRUCT.mat']);
    
            % get position of participant's metadata within metadata
            % variable
            tmp = str2num(erase(subcode{1},'sub-'));
            subs = [];
            for i = 1:size(metadata,2)
                subs = [subs metadata(i).subject];
            end
            subidx = find(subs==tmp)

            % 1. Logistic LASSO

            % configure outputs
            finalcoefs = zeros(size(X,2),1);
            permcoefs = zeros(size(X,2),100);

            % set Y
            Y = [zeros(50,1);ones(50,1)];
            
            % set options structure for training. This includes setting the alpha
            % parameter to 1 (LASSO; default) and lambda to the default range
            % used by WISC_MVPA for SOSLASSO.
            options = glmnetSet;
            options.lambda = linspace(3,0.2);
    
            % train the model on all data
            m = cvglmnet(X,Y,'binomial',options,'deviance',10);

            % extract coefficients
            finalcoefs = m.glmnet_fit.beta(:,find(m.lambda==m.lambda_min));

            for perm = 1:100
                % scramble the data
                Xp = X(PERMUTATION_STRUCT(subidx).permutation_index(:,perm),:);

                % train perm model 
                mp = cvglmnet(Xp,Y,'binomial',options,'deviance',10);
                % extract coefficients using the value of lambda selected by the
                % REAL model
                permcoefs(:,perm) = mp.glmnet_fit.beta(:,find(mp.lambda==m.lambda_min));
            end

            % save results
            save([root,'/work/',subcode{s},'/coefficients/mat/',data{d},'/log.mat'],'finalcoefs','permcoefs');

            % 2. Linear LASSO - 3 dimensions

            for dim = 1:3
                % configure outputs
                finalcoefs = zeros(size(X,2),1);
                permcoefs = zeros(size(X,2),100);
                
                % set options structure for training. This includes setting the alpha
                % parameter to 1 (LASSO; default) and lambda to the default range
                % used by WISC_MVPA for RSL.
                options = glmnetSet;
                options.lambda = linspace(6,0);
        
                % train the model on all data
                m = cvglmnet(X,U(:,dim),'gaussian',options,'deviance',10);
    
                % extract coefficients
                finalcoefs = m.glmnet_fit.beta(:,find(m.lambda==m.lambda_min));
    
                for perm = 1:100
                    % scramble the data
                    Xp = X(PERMUTATION_STRUCT(subidx).permutation_index(:,perm),:);
    
                    % train perm model 
                    mp = cvglmnet(Xp,U(:,dim),'gaussian',options,'deviance',10);
                    % extract coefficients using the value of lambda selected by the
                    % REAL model
                    permcoefs(:,perm) = mp.glmnet_fit.beta(:,find(mp.lambda==m.lambda_min));
                end
    
                % save results
                save([root,'/work/',subcode{s},'/coefficients/mat/',data{d},'/dimension_',num2str(dim),'_lin.mat'],'finalcoefs','permcoefs');
            end
        end
    end
end