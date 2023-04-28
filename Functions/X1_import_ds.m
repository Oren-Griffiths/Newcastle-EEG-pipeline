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
            
%             % find latencies of first and last events.
%             first_times = NaN(1,length(DataConfig.RelevantCodes));
%             last_times = NaN(1,length(DataConfig.RelevantCodes));
% 
%             % find list of all events and latencies
%             allEventTypes = {EEG.event.type}';
%             allEventTimes = {EEG.event.latency}';
% 
%             for thisCode = 1:length(first_times)
%                 % allow for both numeric and string triggers
%                 if isnumeric(DataConfig.RelevantCodes{thisCode})
%                     codeVal = num2str(DataConfig.RelevantCodes{thisCode});
%                     RelEvents_num(thisCode) = DataConfig.RelevantCodes{thisCode};
%                 else
%                     codeVal = DataConfig.RelevantCodes{thisCode};
%                     RelEvents_num(thisCode) = str2double(DataConfig.RelevantCodes{thisCode});
%                 end
%                 strFindResult = strfind(allEventTypes, codeVal);
%                 probeIdx = find(~cellfun(@isempty, strFindResult));
%                 if isempty(probeIdx)
%                 else
%                     first_times(thisCode) = min(probeIdx);
%                     last_times(thisCode) = max(probeIdx);
%                 end
%                 clear probeIdx;
%             end
%             % find the very first instance across events.
%             if ~isnan(min(first_times))
%                 % add a buffer for epoch length before first event
%                 % (EpochMin is neg value, so + fine here).
%                 first_event = allEventTimes{min(first_times)} + 2*DataConfig.EpochMin{1};
%             else
%                 first_event = []; % the very first ms of recording, so no trim.
%             end
%             % find the very last instances across events.
%             if ~isnan(max(last_times))
%                 % add a buffer for epoch length after first event
%                 last_event = allEventTimes{max(last_times)} + 2*DataConfig.EpochMax{1};
%             else
%                 last_event = []; % very last ms. So no trim.
%             end
%             % one more sanity filter...
%             start_time = min(EEG.times);
%             if ~isempty(first_event) && first_event > min(EEG.times)
%                 start_time = first_event;
%             end
%             end_time = max(EEG.times);
%             if ~isempty(last_event) && last_event < max(EEG.times)
%                 end_time = last_event;
%             end
%             % do the start/end trim.
%             EEG = pop_select(EEG, 'time', [start_time/1000, end_time/1000]);

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

%             for thisCode = 1:length(DataConfig.RelevantCodes)
%                 % allow for both numeric and string triggers
%                 if isnumeric(DataConfig.RelevantCodes{thisCode})
%                     codeVal = num2str(DataConfig.RelevantCodes{thisCode});
%                     RelEvents_num(thisCode) = DataConfig.RelevantCodes{thisCode};
%                 else
%                     codeVal = DataConfig.RelevantCodes{thisCode};
%                     RelEvents_num(thisCode) = str2double(DataConfig.RelevantCodes{thisCode});
%                 end
%             end

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
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',[SUB{i} '_ds'],...
            'savenew',[Subject_Path SUB{i} '_ds.set'] ,'gui','off');
        
        % For PREP, only calculate the HEOG and VEOG (no mastoid/average reference reference)
        % choose relevant channel montage.
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
        
        % Add channel location information corresponding to the 3-D coordinates of the electrodes based on 10-10 International System site locations
        EEG = pop_chanedit(EEG, 'lookup',[Current_File_Path filesep 'SupportingDocs' filesep DataConfig.ChanLocs{1}]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_ds_addChans'], ...
            'savenew', [Subject_Path SUB{i} '_ds_addChans.set'], 'gui', 'off');

    end % of subject by subject loop

end