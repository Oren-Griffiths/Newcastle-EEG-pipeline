% Oren's preprocessing script, based on ERP CORE MMN (Luck, Kappenman) and
% core functions of the PrepPipeline toolbox.

% update 15.09.21 to correct some calls to DataConfig.PREP{1} and also to
% make sure ERG channels untouched by cleaning.

function Y1_preprocess_wPREP(DataConfig, SUB)

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
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    %Define subject path based on study directory and subject ID of current subject
    Subject_Path = [DIR filesep SUB{i} filesep];

    % Load the imported and downsampled data.
    FileToOpen = [SUB{i} DataConfig.LastSuffix{1}];
    EEG = pop_loadset( 'filename', FileToOpen, 'filepath', Subject_Path);

    % apply a reference.
    % currently have an average reference. Lets add mastoid reference.
    % Luckily mastoids were passed through cleanline/asr like
    % everything else.
    if strcmp(DataConfig.ReReference{1}, 'Mastoid')
        EEG = pop_reref( EEG, [DataConfig.KeyChans{3},DataConfig.KeyChans{4}] );
    else % average reference
        EEG = pop_reref( EEG, ...
            [DataConfig.firstScalp:DataConfig.lastScalp], ...
            'keepref', 'on');
    end


    % remove any non-scalp channels. With ASR channel removal, hard to
    % keep track of non-scalp channels (which could be removed)
    EEG = pop_select( EEG, ...
        'channel', [DataConfig.firstScalp:DataConfig.lastScalp]);


    % save a copy of the raw data before preprocessing (for later figure).
    rawEEG = EEG;

    %% and now cleanline to remove line noise
    % [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[DataConfig.firstScalp:DataConfig.lastScalp], ...
        'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,...
        'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,...
        'verb',1,'winsize',4,'winstep',1);

    newFilename = [Subject_Path SUB{i} '_ds_cleanline.set'];
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'savenew',newFilename,'gui','off');

    %% and now ASR processing
    % do a highpass filter (to help with later ICA processes)
    EEG = eeg_checkset( EEG );
    % quickly save the details of all the channels, as some are likely
    % about to be removed.
    EEG.urchanlocs = EEG.chanlocs;

    % just remove channels first
    % remove flatlined,  poor correlation or high noise chans
    % do ASR cleaning, but don't reject (just clean)
    % and even then, don't remove periods (via "extra" criterion, box 4)
    % just correct.
    % after that correction is done, itnerpolate the missing channels.
    % and save a file in 'otherData' identifying which were interpreted.

    [EEG,HP,BUR,removed_channels] = clean_artifacts(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4, ...
        'Highpass',[0.25 0.75] ,'BurstCriterion','off','WindowCriterion','off', ...
        'BurstRejection','off','Distance','Euclidian','channels',[] );

    % now we know what channels are missing (removed_channels)
    % finish the ASR algorithm, with as high rank data as possible.
    % (interpolate missing channels later, and reduce rank then).

    % then subject total data set to ASR so as to retain the same number
    % of channels.

    % but need to decide if we're rejecting dirty data, or just cleaning.
    burstRej = 'on';
    if strcmp(DataConfig.ASR_mode, 'interpolate') || strcmp(DataConfig.ASR_mode, 'reject_chan')
        burstRej = 'off';
    end

    EEG = clean_artifacts(EEG, ...
        'FlatlineCriterion','off','ChannelCriterion','off', ...
        'LineNoiseCriterion','off','Highpass','off','BurstCriterion',20, ...
        'WindowCriterion',0.25,'BurstRejection',burstRej,'Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );

    % stash that info away in case we need to reconstruct full data set
    % later.
    if sum(removed_channels) > 0
        % ASR removed some channels. Document which ones and their data.
        EEG.removedChannels.chanlocs = EEG.urchanlocs(removed_channels);
        EEG.removedChannels.data = rawEEG.data(removed_channels,:);

        % done a few things, so check for consistency.
        EEG = eeg_checkset( EEG );
    else % no channels removed. But may need those fields.
        EEG = eeg_checkset( EEG );
        EEG.removedChannels.chanlocs = struct([]); % empty struct
        EEG.removedChannels.data = []; % no data. empty matrix.
    end

    if ~exist([Subject_Path filesep 'Figures' ])
        % create the relevant output subfolders
        mkdir([Subject_Path filesep 'Figures' ]); % Folder for figures
    end
    if ~exist([Subject_Path filesep 'OtherData' ])
        mkdir([Subject_Path filesep 'OtherData' ]); % Folder for summary stats
    end

    % separately save the removed data.
    remChans.chanlocs = EEG.removedChannels.chanlocs;
    remChans.data = EEG.removedChannels.data;
    remChanName = [Subject_Path filesep 'OtherData' filesep SUB{i} '_removedChannels.mat'];
    save(remChanName, "remChans", '-mat');
    clear remChans remChanName;

    %  May need basic SVT applied as 'windowCriterion' in ASR turned off.

    %% only a low pass filter. Not really needed (high pass already done).
    EEG  = pop_basicfilter( EEG,  1:size(EEG.data,1) , 'Boundary', 'boundary', ...
        'Cutoff',  [DataConfig.LPfilter{1}], ...
        'Design', 'butter', 'Filter', 'lowpass', 'Order',  DataConfig.FiltOrder{1}, 'RemoveDC', 'off' );

    %% save the output.
    EEG = pop_saveset( EEG, 'filename', [SUB{i} '_ds_addChans_cleanline_asr_lp_refs.set'],  'filepath', Subject_Path);

    %% and now run some checks and visualize.
    % do simple FFT of pre-filtered data and post-filtered data.
    % condensed steps from mike cohen's examples (with a grid imposed)
    nfft = ceil( EEG.srate/.1 ); % .1 Hz resolution
    hz = linspace(0,EEG.srate,nfft);

    % plot mean FFT across all channels for unfiltered data (rawEEG)
    unfiltFFT = mean(abs( fft(rawEEG.data,nfft,2)/size(rawEEG.data,2)).^2,1);
    filtFFT = mean(abs( fft(EEG.data,nfft,2)/size(EEG.data,2)).^2,1);


    % find index for 50Hz line noise.
    lineNoise = 50;
    [~, idx] = min(abs(hz-lineNoise));
    tolerance = 5; % pull data from 10 bins either side.
    noiseRange = [idx - tolerance, idx + tolerance];
    linenoise_hz = [hz(noiseRange(1)),hz(noiseRange(2))];
    compRange = [ idx - 8*tolerance     , idx - 6*tolerance   ];
    compRange_hz = [hz(compRange(1)),hz(compRange(2))];
    rawRatio = mean(unfiltFFT(noiseRange(1):noiseRange(2)),'all') / ...
        mean(unfiltFFT(compRange(1):compRange(2)),'all');
    filtRatio = mean(filtFFT(noiseRange(1):noiseRange(2)),'all') / ...
        mean(filtFFT(compRange(1):compRange(2)),'all');

    % if lineNoise still 10x larger than relevant comprison range, apply a
    % notch filter directly.
    if filtRatio > 10
        disp('Data still dirty. Applying a notch filter directly.');

        [EEG, com, b] = pop_eegfiltnew(EEG,...
            'locutoff', 49 , ...
            'hicutoff', 51 , ...
            'revfilt',   1 );

        filtFFT_notch = mean(abs( fft(EEG.data,nfft,2)/size(EEG.data,2)).^2,1);
    else
        disp('No additional notch filter needed. Cleanline FTW.');
    end

    % subset the values for plotting
    plot_index = hz<60;
    % plot the values.
    figure;
    line(hz(plot_index), unfiltFFT(plot_index), 'Color', 'black');
    line(hz(plot_index), filtFFT(plot_index),  'Color', 'red');
    if filtRatio > 10
        % and thus we needed a second notch filter.
        line(hz(plot_index), filtFFT_notch(plot_index),  'Color', 'green');
    end
    xline(linenoise_hz(1)   ,  'Color', 'blue', 'LineStyle', '-');
    xline(linenoise_hz(2)   ,  'Color', 'blue', 'LineStyle', '-');
    xline(compRange_hz(1)   ,  'Color', 'blue', 'LineStyle', ':');
    xline(compRange_hz(2)   ,  'Color', 'blue', 'LineStyle', ':');

    % filtered plotted on top.
    % apply some sensible limits.
    y_max = max([ quantile(filtFFT,0.99),  quantile(unfiltFFT,0.99)]);
    if y_max > 0
        ylim([0 y_max]);
    else
        ylim([0 20]);
    end
    title('BandPass. Raw=B,Filt=R.');

    save2pdf([Subject_Path filesep 'Figures' filesep 'X1_' SUB{i} '_BPFilter.pdf']);
    close(gcf);

end % of subject by subject loop

end