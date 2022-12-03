function X5_binEpochs(DataConfig, SUB)

% initialize
close all;
clearvars -except DataConfig SUB;

    % Location of the main study directory
    DIR = fileparts(pwd);
    
    % location of preprocessing files.
    Current_File_Path = pwd;
    
    %**********************************************************************************************************************************************************************
    
    %Loop through each subject listed in SUB
    for i = 1:length(SUB)
        
        %Open EEGLAB and ERPLAB Toolboxes
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        
        %Define subject path based on study directory and subject ID of current subject
        Subject_Path = [DIR filesep SUB{i} filesep];
        
        % Load the continuous ICA-corrected EEG data file 
        FileToOpen = [SUB{i} DataConfig.LastSuffix{1}];
        EEG = pop_loadset( 'filename', FileToOpen, 'filepath', Subject_Path);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_ds_addChans_cleanline_asr_lp_refs_event_weighted'], 'gui', 'off');
        
        % this is now performed in step 1b (preICA), but is retained here
        % in case there's ever a need to bin without precleaning. 
        if size(DataConfig.CustomEpochs) > 0 % custom text epochs to create.
            EventLabels = DataConfig.CustomEpochs;
            DataConfig.CustomEventTriggers = {}; % initialize output cell array.
            for ThisLabel = 1:numel(EventLabels)
                % create a token numeric label per text label. Start at 1001.
                DataConfig.CustomEventTriggers{ThisLabel} = 1000+ThisLabel;
                % cycle through labels and find their appearance
                for ThisEvent = 1:numel(EEG.event)
                    if strcmp (EEG.event(ThisEvent).type,EventLabels(ThisLabel))
                        EEG.event(ThisEvent).type = num2str(1000+ThisLabel);
                    end
                end %cycling vertically through candidate events.
            end % label by label loop
        end
        
         
        %Create EEG Event List containing a record of all event codes and their timing
        EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' }, 'Eventlist', [Subject_Path SUB{i} '_Eventlist.txt'] );
        EEG = eeg_checkset(EEG);
        % must be using numeric codes with Binlister file.
        %Assign events to bins with Binlister; an individual trial may be assigned to more than one bin (bin assignments can be reviewed in each subject's '_Eventlist_Bins.txt' file)
        EEG  = pop_binlister( EEG , 'BDF', [Current_File_Path filesep 'SupportingDocs' filesep  DataConfig.BinListing{1}], 'ExportEL', [Subject_Path SUB{i} '_Eventlist_Bins.txt'], 'IndexEL',  1, 'SendEL2', 'EEG&Text', 'UpdateEEG', 'on', 'Voutput', 'EEG' );
        EEG = eeg_checkset(EEG);
        %Epoch the EEG into 1-second segments time-locked to the response (from -200 ms to 800 ms) and perform baseline correction using the average activity from -200 ms to 0 ms
        EEG = pop_epochbin( EEG , [DataConfig.EpochMin{1} DataConfig.EpochMax{1} ],  [ DataConfig.BaselineMin{1}, DataConfig.BaselineMax{1}]);
        EEG = eeg_checkset(EEG);
        
        %% one more check, which is just for Natalie's reanalysis.
        % one event marker is '32' which is also on the test pin sequence.
        % if that happens, there will be 8 epochs, not 7.
        % the bad epoch will always be first.
        
        % find which epochs contain which markers
        epochvect = cell(1,length(EEG.trials));
        allConds = {'32', 'Condition 32'};
        for k=1:EEG.trials
            [~,t0] = min(abs(cell2mat(EEG.epoch(k).eventlatency)));
            % epochvect(i) = EEG.epoch(i).eventtype{t0};
            epochvect{k} = EEG.epoch(k).eventtype{t0};
        end
        for thisCond = 1:length(allConds)
            if isempty(find(contains(epochvect,allConds{thisCond})))
                index_allConds{thisCond} = 0;
            else
                index_allConds{thisCond} = find(contains(epochvect,allConds{thisCond}));
            end
        end
        % combine all the found epochs
        index = index_allConds{1};
        if length(index_allConds) > 1
            for thisEpoch = 2:length(index_allConds)
                index = union(index, index_allConds{thisEpoch});
            end
        end
        % remove any empty searches
        index = index(index>0);
        % if there's more than one '32' entry, remove the first one as it
        % will be the test pin. 
        if length(index) > 1
            EEG = pop_select( EEG, 'notrial',index(1));
        end
        
        % check to see we haven't broken the format
        EEG = eeg_checkset(EEG);
        % then back to regularly scheduled programming.
        
        %% save and close
        EEG = pop_saveset( EEG, 'filename', [SUB{i} '_ds_addChans_cleanline_asr_lp_refs_event_weighted_epoch_bl.set'], 'filepath', Subject_Path );
        close all;
        
    end % End subject loop
    
end % end of function.
%**********************************************************************************************************************************************************************
