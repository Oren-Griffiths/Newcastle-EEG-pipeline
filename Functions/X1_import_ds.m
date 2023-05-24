% Oren's preprocessing script, based on ERP CORE MMN (Luck, Kappenman) and
% core functions of the PrepPipeline toolbox.

% update 15.09.21 to correct some calls to DataConfig.PREP{1} and also to
% make sure ERG channels untouched by cleaning.

function X1_import_ds(DataConfig, SUB)

% initialize

close all;
clearvars -except DataConfig SUB;

% find the directory one up from this file (for subject folders)
DIR = fileparts(pwd);

% find the directory that this file lives in.
Current_File_Path = pwd;


%% Key input and preprocess loop

%Loop through each subject listed in SUB
for i = 1:length(SUB)

    %Open EEGLAB and ERPLAB Toolboxes
    %Acts as initialization of the relevant variables per participant.
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    %Define subject path based on study directory and subject ID of current subject
    Subject_Path = [DIR filesep num2str(SUB{i}) filesep];

    % record where we're up to in case of crash.
    DataConfig.CurrentSUB = SUB{i};

    % this will almost always be first process called, but just in case
    % it's not let's build in the check to figure out what file to open.

    if isempty(DataConfig.LastProcess)
        % this is the first process for this file, so open the raw data.
        if exist([Subject_Path  SUB{i} '_ds.set']) == 2
            % you've already imported this data set so just start
            % there.
            EEG = pop_loadset( 'filename',[SUB{i} '_ds.set'], 'filepath', Subject_Path);
        else
            switch DataConfig.RawFileType{1}
                case '.bdf'
                    FileToOpen = [Subject_Path  SUB{i} '.bdf'];
                    if isfile(FileToOpen)
                        EEG = pop_biosig(FileToOpen, 'ref', DataConfig.KeyChans{3});
                        EEG = eeg_checkset( EEG );
                    else
                        disp('No .bdf file found');
                        return
                    end
                case '.xdf'
                    FileToOpen = [Subject_Path  SUB{i} '.xdf'];
                    if isfile(FileToOpen)
                        EEG = pop_loadxdf(FileToOpen);
                        EEG = eeg_checkset( EEG );
                    else
                        disp('No .xdf file found');
                        return
                    end
            end
        end

    else
        % not the first file to open. So grab the last file that we
        % processed.
        FileToOpen = [Subject_Path  SUB{i} '.set'];
        if isfile(FileToOpen)
            % shorten FileToOpen to relative address for EEGlab commands
            FileToOpen = [SUB{i} DataConfig.LastSuffix{1}];
            EEG = pop_loadset( 'filename', FileToOpen, 'filepath', Subject_Path);
        else
            MessageForUser = ['No relevant ' DataConfig.LastSuffix{1} ' file found'];
            disp(MessageForUser );
            return
        end
    end

    % if the user supplies "relevant codes" then you can trim and clean
    % the data to only those codes now. Any cleaning here will help both ASR and ICA.
    % if not "RelevantCodes" are supplied, then don't trim.
    if isempty(DataConfig.RelevantCodes)
        % do nothing
    else % start trimming

        % sometimes EEGlab imports numbers (e.g. '23') as simple strings, and
        % sometimes as longer strings (.e.g 'condition 23'). write some code to
        % streamline 'conditionXX' to just 'XX'.
        if ~isempty(DataConfig.RelevantCodes)
            % load a holder variable.
            for thisLabel = 1:length(DataConfig.RelevantCodes)
                % populate cell(strings) and matrix(doubles) for later
                RelCodes_str{thisLabel} = num2str(DataConfig.RelevantCodes{thisLabel});
                RelCodes_num(thisLabel) = DataConfig.RelevantCodes{thisLabel};
                % desired value.
                corrLabel = RelCodes_str{thisLabel} ;
                % add "Condition " (inc space) after each entry.
                wrongLabel = ['condition ' corrLabel ];
                % now replace "condition 23" with "23"
                for ThisEvent = 1:numel(EEG.event)
                    if strcmp (EEG.event(ThisEvent).type,wrongLabel)
                        EEG.event(ThisEvent).type = corrLabel;
                    end
                end %cycling vertically through candidate events.
            end % cycling through relevant codes.
        end % of check

        % noting that event based trimming won't work with xdf string triggers
        % unless the conversion to numeric triggers is done earlier (currently at
        % step 3).


        % and do an event-based trim within that set too.
        EEG = pop_erplabDeleteTimeSegments(EEG, ...
            'timeThresholdMS', 2 * max(abs(DataConfig.EpochMin{1}), abs(DataConfig.EpochMax{1}) ), ...
            'beforeEventcodeBufferMS', 2*abs(DataConfig.EpochMin{1}), ...
            'afterEventcodeBufferMS', 2*abs(DataConfig.EpochMax{1}), ...
            'ignoreUseEventcodes', RelCodes_num, ...
            'ignoreUseType', 'use',  ...
            'ignoreBoundary', 1  );

    end % event-based trimming.

    if DataConfig.DownSample{1} == EEG.srate
        % already at the desired sample rate so leave it be.
    else
        % Downsample from the recorded sampling rate e.g. 2048 Hz to e.g. 256 Hz
        % to speed data processing (automatically applies the appropriate
        % low-pass anti-aliasing filter)
        EEG = pop_resample(EEG, DataConfig.DownSample{1});
    end
    
    % save downsampled, or non-downsampled, as new data set.
    pop_saveset(EEG, 'filename', [SUB{i} '_ds.set'] , 'filepath', Subject_Path);


    % if no information came from the import, then add channel labels.
    if isempty(EEG.chanlocs)
        switch DataConfig.RawFileType{1}
            case '.bdf'
                NoRefChanLbls = [Current_File_Path filesep 'SupportingDocs' ...
                    filesep 'ChannelsFor' num2str(DataConfig.TotalChannels{1}) '_NoRef_BDF.txt'];
            case '.xdf'
                NoRefChanLbls = [Current_File_Path filesep 'SupportingDocs' ...
                    filesep 'ChannelsFor' num2str(DataConfig.TotalChannels{1}) '_NoRef_XDF.txt'];
        end
        % apply that montage.
        EEG = pop_eegchanoperator( EEG, NoRefChanLbls, 'Saveas', 'off');
        disp('Applied manual channel labels.')
    end
    
    % Add channel location information corresponding to the 3-D coordinates of the electrodes based on 10-10 International System site locations
    EEG = pop_chanedit(EEG, 'lookup',[Current_File_Path filesep 'SupportingDocs' filesep DataConfig.ChanLocs{1}]);
    pop_saveset(EEG, 'filename', [SUB{i} '_ds_addChans.set'] , 'filepath', Subject_Path);

end % of subject by subject loop

end