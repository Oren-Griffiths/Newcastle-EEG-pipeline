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
        
        % earlier call before more flexible opening command included.
        % FileToOpen = [Subject_Path  SUB{i} '.bdf'];
        
        % Take the .bdf file and make it a .set file (Oren added this)
        % Token reference is channel 40 for input, and then will reference
        % properly in a moment. Files must be named: "23.bdf" on the way in
        % and will be named "17.set" on the way out.
        
        % If file are e.g. "P17.bdf" then just adjust SUB cell array with
        % appropriate values.
        
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