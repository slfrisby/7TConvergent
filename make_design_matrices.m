% by Saskia. To make design matrices for decoding, specific to the order in
% which each participant saw the stimuli.

mkdir('/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/')
cd('/imaging/projects/cbu/wbic-p00591-DAISY/main/behavioural/')
load('/imaging/projects/cbu/wbic-p00591-DAISY/main/scripts/stimulimaster.mat')

% initialise variables
badsubj = [];

% loop over participants
for s = 1:32
    
    % set subcode
    subcode = sprintf('sub-%03d',s);
    
    % check that there are 4 csvs of stimuli, one for each run. If there
    % are missing runs, don't make design matrices for that participant (but store
    % the subcodes in a string array)
    
    % list all files in csv directory beginning with the subcode
    tmp = dir([pwd,'/csv/',sprintf('%03d',s),'*']);
    
    % if there are missing files
    if size(tmp,1) ~= 4
        % store the subcode
        badsubj = [badsubj,convertCharsToStrings(subcode)];
        % move on to the next participant
    % if there are no missing files
    elseif size(tmp,1) == 4
        % go onto the next part of the analysis.
        
        % load table of accuracies
        excelfile = readtable('accuracies.xlsx','Sheet',sprintf('%03d',s));
        
        % loop over runs 
        for r = 1:4
            
            % extract the column of the excel file containing stimulus
            % names. (The y-coordinate is a sequence that happens to give
            % the right column in the table!)
            stimuli = excelfile(1:100,3*r-2);
            
            % initialise components of design matrix
            names = cell(1,100);
            onsets = cell(1,100);
            durations = cell(1,100);
            
            % loop over stimuli
            % get rid of unwanted text in stimuli names
            for stim = 1:size(stimuli,1)
                stimname = char(stimuli{stim,1});
                a = erase(stimname,'StimFiles/');
                b = erase(a,'.bmp');
                names{stim} = b;
                
                % set durations - 1.5 seconds per stimulus (based on
                % reaction time norms)
                durations{stim} = 1.5;
                
                % set onsets 
                onsets{stim} = 16 + 8*(stim-1);
            end
            % save design matrix for decoding
            save(['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode,'_run-0',num2str(r),'_type-multivariate_design-matrix.mat'],'names','onsets','durations')
            
            % create design matrices for univariate. First, stimulus vs. fixation. 
            % initialise
            animatecorrect = [];
            inanimatecorrect = [];
            animateincorrect = [];
            inanimateincorrect = [];
            
            % loop over stimuli
            for stim = 1:size(stimuli,1)
                % if the stimulus is animate (i.e. it is in the first half
                % of the stimuli master)
                if find(strcmp(names{stim},stimulimaster)) <= 50
                    % and the person named it correctly
                    if table2array(excelfile(stim,3*r-1))
                       % add its onset time to the list of animate correct
                       % onsets
                       animatecorrect = [animatecorrect,onsets{stim}];
                    % else if the person did not name it correctly
                    elseif ~table2array(excelfile(stim,3*r-1))
                       % add its onset time to the list of animate
                       % incorrect onsets
                       animateincorrect = [animateincorrect, onsets{stim}];
                    end
                    % else if the stimulus is inanimate
                elseif find(strcmp(names{stim},stimulimaster)) > 50
                    % and the person named it correctly
                    if table2array(excelfile(stim,3*r-1))
                       % add its onset time to the list of inanimate correct
                       % onsets
                       inanimatecorrect = [inanimatecorrect,onsets{stim}];
                    % else if the person did not name it correctly
                    elseif ~table2array(excelfile(stim,3*r-1))
                       % add its onset time to the list of inanimate
                       % incorrect onsets
                       inanimateincorrect = [inanimateincorrect, onsets{stim}];
                    end
                end
            end
          
            % construct design matrix
            names = {'AnimateCorrect','InanimateCorrect','AnimateIncorrect','InanimateIncorrect'};
            durations = {1.5 1.5 1.5 1.5};
            onsets = {animatecorrect,inanimatecorrect,animateincorrect,inanimateincorrect};

            % sometimes a participant does not make any errors. In this
            % case, SPM will throw an error because the condition will have
            % no onsets. If this is the case, set the onset for this
            % condition to be 10000 (a time that does not exist).
            for cond = 1:4
                if isempty(onsets{cond})
                    onsets{cond} = 10000;
                end
            end

            % save
            save(['/imaging/projects/cbu/wbic-p00591-DAISY/main/work/design_matrices/',subcode,'_run-0',num2str(r),'_type-univariate_design-matrix.mat'],'names','onsets','durations')
        end
    end
end