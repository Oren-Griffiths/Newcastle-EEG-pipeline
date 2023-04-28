
function X4_AutoRunTheICA(DataConfig,SUB, ICAmode)

% initialize
close all;
clearvars -except DataConfig SUB ICAmode;

    % Location of the main study directory
    DIR = fileparts(pwd);
    
    % location of preprocessing files.
    Current_File_Path = pwd;
    
    %***********************************************************************************************************************************************
    
    %Loop through each subject listed in SUB
    for i = 1:length(SUB)
        
        %Open EEGLAB and ERPLAB Toolboxes
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        
        %Define subject path based on study directory and subject ID of current subject
        Subject_Path = [DIR filesep SUB{i} filesep];
        
        % record where we're up to in case of crash.
        DataConfig.CurrentSUB = SUB{i};
        
        % Load the semi-continuous EEG data file outputted from Script #2 in .set EEGLAB file format
        FileToOpen = [SUB{i} DataConfig.LastSuffix{1}];
        EEG = pop_loadset( 'filename', FileToOpen, 'filepath', Subject_Path);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_preICA'], 'gui', 'off');
        
        % runica can sometimes miscalculate rank and screw up. So make it
        % guess the rank, and force that rank onto the calculation using
        % 'pca' key value pair to avoid complex components that break things.
        ChansForICA = 1:size(EEG.data, 1); % forcing ICA to all channels. -
        EffectiveRank = rank(EEG.data(ChansForICA,:));
        
        EEG = pop_runica(EEG,'extended',0,'chanind', ChansForICA, 'pca',EffectiveRank);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_ds_addChans_cleanline_asr_lp_refs_event_weighted'], 'savenew', [Subject_Path SUB{i} '_ds_addChans_cleanline_asr_lp_refs_event_weighted.set'], 'gui', 'off');
        
        % new, an automated ICA removal using ICLabel plug in.
        % Using ICLabel output is more about removing components.
        % So it is placed in  process X4 instead.

        EEG = iclabel(EEG, 'default');

        % save the newly complete file to disk.
        EEG = pop_saveset( EEG, 'filename', [SUB{i} '_ds_addChans_cleanline_asr_lp_refs_event_weighted.set'],'filepath', Subject_Path);

        % Pick the relevant componenets to remove/keep. 
        ic_threshold = 0.5; % 50% makes sense here, as must be most likely classification because
        % metric of ICLabel output is posterior probability (e.g. sum to 1).
        brainIdx  = find(EEG.etc.ic_classification.ICLabel.classifications(:,1) >= ic_threshold);
        eyeIdx = find(EEG.etc.ic_classification.ICLabel.classifications(:,3) >= ic_threshold);
        
        % remove the components and keep a tab of what was removed/kept. 
        if strcmp(ICAmode, 'keepBrain')
            save([Subject_Path 'OtherData' filesep 'BrainComponentsKept.mat'], 'brainIdx');
            EEG = pop_subcomp(EEG, brainIdx, 0, 1); % updates IClabel fields too.
        elseif strcmp(ICAmode, 'removeEyes')
            save([Subject_Path 'OtherData' filesep 'EyeComponentsRejected.mat'], 'eyeIdx');
            EEG = pop_subcomp(EEG, eyeIdx, 0, 0); % updates IClabel fields too.
        else
            disp('Misspecified ICA behaviour. No components removed');
        end
        
        % draw and save topoplots per components.
        figure;
        pop_topoplot(EEG, 0, [1:size(EEG.icaweights,1)],[SUB{i} '_ds_addChans_cleanline_asr_lp_refs_event_weighted'], [round(sqrt(size(EEG.icaweights,1)))+1 round(sqrt(size(EEG.icaweights,1)))+1] ,0,'electrodes','on');
        % add some info about which components the algorith wants out/in.
        % (if available). 
        if isfield(EEG.etc, 'ic_classification')
            eyeIdx = find(EEG.etc.ic_classification.ICLabel.classifications(:,3) >= 0.5);
            componentsText = ['X3b_Eyes_' num2str(eyeIdx') '_SUB_' SUB{i} '_ICA_Weights.pdf'];
        else
            componentsText =  ['X3b_SUB_' SUB{i} '_ICA_Weights.pdf'];
        end
        save2pdf([Subject_Path filesep 'Figures' filesep componentsText]);
        close all

        % interpolate any missing channels from ASR process.
        if strcmp(DataConfig.ASR_mode, 'interpolate') || strcmp(DataConfig.ASR_mode, 'reject_time')

            % load in the removed channels.
            remChanFile = [Subject_Path filesep 'OtherData' filesep SUB{1} '_removedChannels.mat'];
            if exist(remChanFile) > 0
                load(remChanFile);
                % loads a struct variable called remChans
                % with fields chanlocs and data
                % empty if no chans were removed.
                if ~isempty(remChans.chanlocs)
                    % ok, let's do some interpolation.
                    numChanRem = numel(remChans.chanlocs); % number removed
                    numChanRetained = size(EEG.data, 1); % number retained
                    disp('Interpolating ASR removed channels.');
                    disp(['Retained channels:' num2str(numChanRetained) '.']);
                    disp(['Removed channels:' num2str(numChanRem) '.']);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    placeholder_data = ones(numChanRem ,size(EEG.data,2));
                    EEG.chanlocs = [EEG.chanlocs, remChans.chanlocs];
                    EEG.data = [EEG.data; placeholder_data];
                    EEG.nbchan = EEG.nbchan + numChanRem;
                    % do the spherical interpolation
                    EEG = pop_interp(EEG, [EEG.nbchan - (numChanRem-1): EEG.nbchan], 'spherical');
                    % ok, we'll need to reshape data back into original order.
                    % for next steps, move from struct into table
                    t_chanloc = struct2table(EEG.chanlocs);
                    data_index = 1:height(t_chanloc);
                    t_chanloc.data_index = [data_index'];
                    % reorder the chanloc array
                    t_chanloc = sortrows(t_chanloc,'urchan');
                    % export the reordering index
                    data_index = t_chanloc.data_index;
                    % now delete the extra chanloc 0field, and revert to structure.
                    t_chanloc = removevars(t_chanloc, 'data_index');
                    EEG.chanlocs = table2struct(t_chanloc);

                    % now reorder the raw data too.
                    for thisChan = 1:size(EEG.data,1)
                        EEG.data(thisChan,:) =  EEG.data(data_index(thisChan),:);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    disp('Found removed channels file, but no channels removed');
                end
            else
                disp('Cannot find removed channels file in OtherData');
            end % checking for removed channels file
        end % of interpolation mode check.

        % save the post_ICA and interpolation file to disk.
        EEG = pop_saveset( EEG, 'filename', [SUB{i} '_ds_addChans_cleanline_asr_lp_refs_event_weighted.set'],'filepath', Subject_Path);

    end % End subject loop

end
%***********************************************************************************************************************************************
