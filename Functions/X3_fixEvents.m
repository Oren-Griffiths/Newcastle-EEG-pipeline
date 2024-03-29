
function X3_fixEvents(DataConfig, SUB)
% initialize
close all;
clearvars -except DataConfig SUB;

% Location of the main study directory
DIR = fileparts(pwd);

% location of preprocessing files.
Current_File_Path = pwd;

%% events to correct and how to adjust them.
% some events need to be moved.

% Key input and preprocess loop
%Loop through each subject listed in SUB
for i = 1:length(SUB)
    
    % Open EEGLAB and ERPLAB Toolboxes
    % Acts as initialization of the relevant variables per participant.
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    
    % Define subject path based on study directory and subject ID of current subject
    Subject_Path = [DIR filesep num2str(SUB{i}) filesep];
    
    % record where we're up to in case of crash.
    DataConfig.CurrentSUB = SUB{i};
    
    % open the relevant file and get cracking.
    % should be an earlier file to load, but just in case there's not.
    if isempty(DataConfig.LastProcess)
        % this is the first process for this file, so open the raw data.
        switch DataConfig.RawFileType{1}
            case '.bdf'
                FileToOpen = [SUB{i} '.bdf'];
                EEG = pop_biosig([Subject_Path  FileToOpen], 'ref', DataConfig.KeyChans{5});
                EEG = eeg_checkset( EEG );
            case '.xdf'
                FileToOpen = [SUB{i} '.xdf'];
                EEG = pop_loadxdf([Subject_Path  FileToOpen]);
                EEG = eeg_checkset( EEG );
        end
    else
        % not the first file to open. So grab the last file that we
        % processed.
        FileToOpen = [SUB{i} DataConfig.LastSuffix{1}];
        EEG = pop_loadset( 'filename', FileToOpen, 'filepath', Subject_Path);
    end
    
    %Load the continuous EEG data file outputted from Script #1a in .set EEGLAB file format
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_ds_addChans_PREP_bp_refs'], 'gui', 'off');

    % if triggers are text (e.g. LSL) then adjust them here (so can
    % use "relevant codes" to limit ICA windows).
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
    
    if isempty(DataConfig.EventsToAdjust)
        % then don't adjust any events. Just save with an updated file
        % name.
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_ds_addChans_PREP_bp_refs_event'], 'savenew', [Subject_Path SUB{i} '_ds_addChans_PREP_bp_refs_event.set'], 'gui', 'off');
    else % else shift the specified events by the specified amount.
        % check to see that first nominated event was moved as
        % intended. (Assume the rest follow ok).
        
        % if there are edftype values that match our target set,
        % push them into type field.
        setOfEvents = cell2mat(DataConfig.RelevantCodes);
        for thisEvent = 1:numel(EEG.event)
            if isempty(intersect(EEG.event(thisEvent).edftype, setOfEvents))
                % do nothing for irrelevant values.
            else % if an edftype matches, push that into type field.
                EEG.event(thisEvent).type = num2str(EEG.event(thisEvent).edftype);
                disp(['Found one here: ' num2str(EEG.event(thisEvent).urevent)]);
            end
        end
        
        preEvents = struct2table(EEG.event);
        % sometimes this will change edftype into cell of strings,
        % sometimes into cell of numerics. Need to force string.
        if strcmp(class(preEvents.edftype{1}), 'double')
            % it's made a cell of numbers, so force a change.
            preEvents.edftype = cellfun(@num2str, preEvents.edftype, 'UniformOutput', false);
            preEvents.type = cellfun(@num2str, preEvents.type, 'UniformOutput', false);
        end
        
        % now adjust the specified events.
        % should just use the built in EEGlab function.
        EEG = pop_adjustevents(EEG,'addms', DataConfig.TimeShift{1}*1000, 'eventtypes', DataConfig.EventsToAdjust);
        % grab the events after the switch.
        postEvents = struct2table(EEG.event);
        % sometimes this will change edftype into cell of strings,
        % sometimes into cell of numerics. Need to force string.
        if strcmp(class(postEvents.edftype{1}), 'double')
            % it's made a cell of numbers, so force a change.
            postEvents.edftype = cellfun(@num2str, postEvents.edftype, 'UniformOutput', false);
            postEvents.type = cellfun(@num2str, postEvents.type, 'UniformOutput', false);
        end
        
        for ThisBin = 1:numel(DataConfig.EventsToAdjust)
            % need to do a little extra here becuase sometimes the
            % relevant code goes into edftype field and sometimes
            % it goes into type field. No clear way to know
            % pre-import, so we'll just take both.
            
            % find pre events for ThisBin
            indxtmp1 = strcmp(preEvents.type, DataConfig.EventsToAdjust{ThisBin});
            indxtmp2 = strcmp(preEvents.edftype,DataConfig.EventsToAdjust{ThisBin});
            indxtmp = or(indxtmp1, indxtmp2);
            preEventsPerBin = preEvents(indxtmp,:);
            
            % find post events for ThisBin
            indxtmp1 = strcmp(postEvents.type, DataConfig.EventsToAdjust{ThisBin});
            indxtmp2 = strcmp(postEvents.edftype,DataConfig.EventsToAdjust{ThisBin});
            indxtmp = or(indxtmp1, indxtmp2);
            postEventsPerBin = postEvents(indxtmp,:);
            
            % compare the pre and post events and output the change.
            heightDiff = height(preEventsPerBin) - height(postEventsPerBin);
            if heightDiff == 0
                % then they're the same length, so just compare
                % directly.
                latencyDiff = preEventsPerBin.latency - postEventsPerBin.latency;
                % initially measured in samples, so convert to sec.
                latencyDiff = latencyDiff./EEG.srate;
                figure;
                plot(latencyDiff);
                title(['Adjustment of Bin ' num2str(ThisBin) ' in seconds']);
                saveas(gcf, [Subject_Path filesep 'Figures' filesep 'X1b_' SUB{i} ...
                    '_Bin_' num2str(ThisBin) '_EventsMoved.png']);
                close(gcf);
            else if heightDiff > 0
                    % then one or more events were chopped off pre to
                    % post.
                    preEventsPerBin = preEventsPerBin(heightDiff+1:end, :);
                    latencyDiff = preEventsPerBin.latency - postEventsPerBin.latency;
                    % initially measured in samples, so convert to sec.
                    latencyDiff = latencyDiff./EEG.srate;
                    figure;
                    plot(latencyDiff);
                    title(['Adjustment of Bin ' num2str(ThisBin) ' in seconds']);
                    saveas(gcf, [Subject_Path filesep 'Figures' filesep  'X1b_' SUB{i} ...
                        '_Bin_' num2str(ThisBin) '_EventsMoved.png']);
                    close(gcf);
                else % something has gone wrong. Can't be adding events.
                    disp('Error: Somehow events were added when shifting times.');
                end
            end % of check for different table heights.
        end % of bin by bin loop.
        
    end % of the part to skip if no events to adjust
    % format for creating a new .set file and saving it.
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_ds_addChans_cleanline_asr_lp_refs_event'], 'savenew', [Subject_Path SUB{i} '_ds_addChans_cleanline_asr_lp_refs_event.set'], 'gui', 'off');
    
    % and even if it is skipped, output a table of the number of events
    % in a .txt file. But exclude 0,1, 256 codes which are just noise.
    disp('Calculating and reporting event tallies.');
    eventTable = struct2table(EEG.event);
    % convert to numeric and add a column called "codes"
    eventCodes = str2double(eventTable.type);
    eventTable.codes = eventCodes;
    % filter 256, 0, 1, values.
    nanCodes = isnan(eventCodes); % remove text codes
    badCodes = eventCodes == 1 | eventCodes ==0 | eventCodes ==256;
    % removes pointless 0, 1, 256 auto-generated codes.
    goodCodes = ~(nanCodes | badCodes);
    eventTable = eventTable(goodCodes,:);
    % remove unnecessary fields.
    eventTable = table(eventTable.codes, 'VariableNames', {'codes'});
    if height(eventTable) < 1
        % for text triggers (in .xdf files) there will be no
        % "goodCodes" as all text will be coded as a NaN. So no point
        % doing tallies. Just skip over everything.
    else
        for ThisRow = 1:height(eventTable)
            if eventTable.codes(ThisRow) > 256
                eventTable.corrCodes(ThisRow) = eventTable.codes(ThisRow) - 256;
            else
                eventTable.corrCodes(ThisRow) = eventTable.codes(ThisRow);
            end
        end
        
        % count the unique values and report them as .csv files.
        [C,ia,ic] = unique(eventTable.codes);
        code_counts = accumarray(ic,1);
        countsPerCode = array2table([C, code_counts]);
        countsPerCode.Properties.VariableNames(1:2) = {'event','count'};
        outFileName = [Subject_Path SUB{i} '_eventCounts_initial.csv'];
        writetable(countsPerCode, outFileName);
        % do the same, but with 256 adjustment.
        [C,ia,ic] = unique(eventTable.corrCodes);
        code_counts = accumarray(ic,1);
        countsPerCorrCode = array2table([C, code_counts]);
        countsPerCorrCode.Properties.VariableNames(1:2) = {'event','count'};
        outFileName= [Subject_Path SUB{i} '_eventCountsWith256Correction.csv'];
        writetable(countsPerCorrCode, outFileName);
    end % of skipping over tallies when text codes are used.
end % of subject by subject loop

% record the last process performed. (pointless now with parallel
% execution as DataConfig is not a global or passed out).

end % of function.