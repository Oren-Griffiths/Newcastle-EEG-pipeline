% Oren's preprocessing script, based on ERP CORE MMN (Luck, Kappenman) and
% core functions of the PrepPipeline toolbox.

% update 15.09.21 to correct some calls to DataConfig.PREP{1} and also to
% make sure ERG channels untouched by cleaning.

function Y1_preprocess_wPREP(DataConfig, SUB, mode)

% initialize

close all;
clearvars -except DataConfig SUB mode;

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
       EEG = eeg_checkset( EEG );
       
       % just remove channels first
       [EEG,HP,BUR,removed_channels] = clean_artifacts(EEG, ...
           'FlatlineCriterion',5,'ChannelCriterion',0.8, ...
           'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,'BurstCriterion','off', ...
           'WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
       
       % stash that info away in case we need to reconstruct full data set
       % later.
       if sum(removed_channels) > 0
           % ASR removed some channels. Document which ones and their data.
           EEG.removedChannels.chanlocs = EEG.urchanlocs(removed_channels);
           EEG.removedChannels.data = rawEEG.data(removed_channels,:);
           % (this will only work for this one data set...)
           DataConfig.lastScalp = DataConfig.lastScalp - numel(EEG.removedChannels.chanlocs);
           
           % interpolate, if that's requested.
           if strcmp(mode, 'interpolate')
               EEG.chanlocs = [EEG.chanlocs, EEG.removedChannels.chanlocs];
               EEG.data = [EEG.data; EEG.removedChannels.data];
               EEG.nbchan = EEG.nbchan + numel(EEG.removedChannels.chanlocs);
               % do the spherical interpolation
               EEG = pop_interp(EEG, [size(EEG.data) - (numel(EEG.removedChannels.chanlocs)-1): size(EEG.data,1)], 'spherical');
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
               
           end
           % done a few things, so check for consistency.
           EEG = eeg_checkset( EEG );
       else % no channels removed. But may need those fields. 
           EEG.removedChannels.chanlocs = struct([]); % empty struct
           EEG.removedChannels.data = []; % no data. empty matrix.
       end
       
       % then finish the ASR algorithm, now we have the desired 
       % number of channels.
       
       % then subject total data set to ASR so as to retain the same number
       % of channels. 
       EEG = clean_artifacts(EEG, ...
           'FlatlineCriterion','off','ChannelCriterion','off', ...
           'LineNoiseCriterion','off','Highpass','off','BurstCriterion',20, ...
           'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian', ...
           'WindowCriterionTolerances',[-Inf 7] );

       newFilename = [Subject_Path SUB{i} '_ds_cleanline_asr.set'];
       EEG = eeg_checkset( EEG );
       [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[SUB{i} '_ds_addChans_cleanline_asr'],'savenew',newFilename,'gui','off');
        
        %% only a low pass filter. Not really needed (high pass already done). 
        EEG  = pop_basicfilter( EEG,  1:size(EEG.data,1) , 'Boundary', 'boundary', ...
            'Cutoff',  [DataConfig.LPfilter{1}], ...
            'Design', 'butter', 'Filter', 'lowpass', 'Order',  DataConfig.FiltOrder{1}, 'RemoveDC', 'off' );
        % save at the bandpass point.
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_ds_addChans_cleanline_asr_lp'], ...
            'savenew', [Subject_Path SUB{i} '_ds_addChans_cleanline_asr_lp.set'], 'gui', 'off');

        %% save the output.
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_ds_addChans_cleanline_asr_lp_refs'], ...
            'savenew', [Subject_Path SUB{i} '_ds_addChans_cleanline_asr_lp_refs.set'], 'gui', 'off');
        
        %% and now run some checks and visualize.
        % do simple FFT of pre-filtered data and post-filtered data. 
        % condensed steps from mike cohen's examples (with a grid imposed)
        nfft = ceil( EEG.srate/.1 ); % .1 Hz resolution
        hz = linspace(0,EEG.srate,nfft);
        
        % plot mean FFT across all channels for unfiltered data (rawEEG)
        unfiltFFT = mean(abs( fft(rawEEG.data,nfft,2)/size(rawEEG.data,2)).^2,1);
        filtFFT = mean(abs( fft(EEG.data,nfft,2)/size(EEG.data,2)).^2,1);

        % if you need to do some debugging of drawing, save workspace now.
        % save('Debug_workspace.mat');
        
        if ~exist([Subject_Path filesep 'Figures' ])
            % create the relevant output subfolders
            mkdir([Subject_Path filesep 'Figures' ]); % Folder for figures
            mkdir([Subject_Path filesep 'OtherData' ]); % Folder for summary stats
        end

        % subset the values for plotting
        plot_index = hz<60;
        % plot the values.
        figure;
        line(hz(plot_index), unfiltFFT(plot_index), 'Color', 'k');
        line(hz(plot_index), filtFFT(plot_index),  'Color', 'r');
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