%% input parameters.

% for gadi to access eeglab.
addpath('/home/575/og5989/EEGlab/eeglab2022.1');

allConditions = {'B1(' , 'B2(' , 'B3(' , ...
    'B4(' , 'B5(' , 'B6(', 'B7(', 'B8(', 'B9('};
% do these conditions have different baselines? If so, you'll need to
% complete the section below (1 = custom, 0 = shared).
customBaseline = 1;

% how many threads to use in parallel mode? Default is numThreads = 1.
numThreads = 4;

% which channels do you want to use?
% currently written to loop through all scalp channels.
% keyChans = 1:32; % now it just uses all the available scalp channels.

%% header structure grabs file and config data

% what's the relevant config file called?
ConfigFileName = 'Config_CaitlinCharlotte_v2023';

Current_File_Path = pwd;
addpath('Functions');
ConfigFilePath = [Current_File_Path filesep 'SupportingDocs' filesep ConfigFileName '.xlsx'];
Options = detectImportOptions(ConfigFilePath);

for k = 1:numel(Options.VariableTypes)
    Options.VariableTypes{k} = 'char';
end
DataConfig = table2struct(readtable(ConfigFilePath, Options));
DataConfig = adjustConfigData(DataConfig);

% and open eeglab to access the EEGlab functions
eeglab;
% just shorten variable name
SUB = DataConfig.SUB;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manual override for troubleshooting.

%% add any custom baseline values.%%%%%%%%%%%%%%%%%%%%%

if customBaseline == 1
    % then we need values for each of the conditions declared.
    binTimings = struct;
    % add token values which can then be overwritten by manual values.
    for k = 1:size(allConditions, 2)
        binTimings(k).baseline = [-Inf, 0];
        binTimings(k).measureWindow = [0, Inf];
        binTimings(k).label = allConditions{k};
    end
    % custom values can then be overwritten as required in section below.

    % all these timings are relative to the commencement of the 8700ms
    % maintenance period.
    binTimings(1).baseline = [-6500, -5500];
    binTimings(1).measureWindow = [0 8700];
    binTimings(1).label = 'B1';
    %
    binTimings(2).baseline = [-6500, -5500];
    binTimings(2).measureWindow = [0 8700];
    binTimings(2).label = 'B2';
    %
    binTimings(3).baseline = [-6500, -5500];
    binTimings(3).measureWindow = [0 8700];
    binTimings(3).label = 'B3';
    %
    binTimings(4).baseline = [-6500, -5500];
    binTimings(4).measureWindow = [0 8700];
    binTimings(4).label = 'B4';
    %
    binTimings(5).baseline = [-7535, -6535];
    binTimings(5).measureWindow = [0, 8700];
    binTimings(5).label = 'B5';
    %
    binTimings(6).baseline = [-12282 -11282];
    binTimings(6).measureWindow = [0, 8700];
    binTimings(6).label = 'B6';
    %
    binTimings(7).baseline = [-7535, -6535];
    binTimings(7).measureWindow = [0, 8700];
    binTimings(7).label = 'B7';
    %
    binTimings(8).baseline = [-12282 -11282];
    binTimings(8).measureWindow = [0, 8700];
    binTimings(8).label = 'B8';
    %
    binTimings(9).baseline = [-6000 -5000];
    binTimings(9).measureWindow = [0, 8700];
    binTimings(9).label = 'B9';
end

%% start processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize an output folder, if necessary
if exist('TF_output', 'dir') == 7
else
    mkdir 'TF_output'
end

parfor (k = 1:length(SUB), numThreads)

    PIDfolder = [fileparts(pwd) filesep SUB{k}];
    inputFile = [PIDfolder filesep SUB{k} '_ds_addChans_cleanline_asr_lp_refs_event_weighted_epoch_bl_ar.set'];

    %% find each file and open it up
    EEG = pop_loadset('filename',  inputFile);

    %  if input files are missing channel information for some reason, you can add it now.
    % EEG = pop_chanedit(EEG, 'lookup',[Current_File_Path filesep 'SupportingDocs' filesep DataConfig.ChanLocs{1}]);

    % initialize an output variable for the collection of tf_data in parallel.
    TotalTF = struct;
    TotalTF.PID = SUB{k};
    TotalTF.inputFile = inputFile;
    TotalTF.inputFolder = PIDfolder;
    TotalTF.data = [];

    %% separate  them out so each data set = 1 bin type.

    % first, select the epochs that we want.
    % populate a list with of the time-locking event per epoch.

    for thisCond = 1:length(allConditions)

        cond2use = allConditions{thisCond};

        % report progress.
        disp(['Starting with condition ' cond2use ' in PID ' SUB{k}])

        epochvect = cell(1,EEG.trials);
        keepTrials = zeros(1,EEG.trials);
        for i=1:EEG.trials
            [~,t0] = min(abs(cell2mat(EEG.epoch(i).eventlatency)));
            % epochvect(i) = EEG.epoch(i).eventtype{t0};
            epochvect{i} = EEG.epoch(i).eventtype{t0};
            % if strcmp(epochvect{i},cond2use)
            if contains(epochvect{i},cond2use)
                keepTrials(i) = 1;
            end
        end

        % limit our data to those epochs.
        data  = EEG.data(:,:,logical(keepTrials));

        %%  prepare parameters for repeated tf call.
        frames = EEG.pnts;
        tlimits = [EEG.times(1) EEG.times(end)];
        srate = EEG.srate;
        % cycles = [3 0.5]; % standard parameters for now.
        cycles = [1 0.5]; % allows us to get down into theta range.
        % cycles = 0; % FFT with Hanning tapering (no wavelets)

        % move these to header later.
        keyChans = 1:size(EEG.data,1);
        freqRange = [0 30];
        tempResolution = 200; % how many points to output on x-axis. (44s data, so 440 = 0.1s)
        baseline = []; % baseline period measured in epoch time (ms)

        % time-freq for data, per channel.
        % initialise some varialbes
        tf_data = struct;
        tf_data = []; % to establish order independence for parfor.
        tf_data.times = [];
        tf_data.freqs = [];

        for thisChan = keyChans

            [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = ...
                newtimef(data(thisChan,:,:),EEG.pnts, [tlimits], 256, cycles, ...
                'plotersp', 'off', 'plotitc', 'off', ...
                'baseline', binTimings(thisCond).baseline, ...
                'trialbase', 'on', ...
                'freqs', [1, 40], ...
                'verbose', 'off');

            tf_data.PID = SUB{k};
            tf_data.chan(thisChan).lbl = EEG.chanlocs(thisChan).labels;
            tf_data.chan(thisChan).ersp = ersp;
            tf_data.chan(thisChan).itc = abs(itc);

            % check to see if you need to write before writing.
            if isempty(tf_data.times)
                tf_data.times = times;
            end

            if isempty(tf_data.freqs)
                tf_data.freqs = freqs;
            end

        end % of channel by channel loop

        TotalTF.data.cond(thisCond).chan = tf_data.chan;
        TotalTF.data.cond(thisCond).times = tf_data.times;
        TotalTF.data.cond(thisCond).freqs = tf_data.freqs;
        % report progress.
        disp(['Finished with condition ' cond2use ' in PID ' SUB{k}])

    end % of condition by condition loop

outName =  [pwd filesep 'TF_output' filesep SUB{k} '_TFdata.mat'];
disp(['Saving file: ' outName])
parsave(outName, TotalTF );

end % of PID looping cycle
