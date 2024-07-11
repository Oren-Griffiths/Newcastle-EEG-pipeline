
% parameters
% keyChans = {'FC1', 'FC3', 'FC5', 'C1', 'C3', 'C5', 'F1', 'F3', 'F5'};
keyChans = {'all'};
% use 'all' for all scalp channels (i.e. global power).

%keyChans = [5, 26, 31, 4, 27, 2, 29];
% for a set of 7 sensors centred on Fz (32 chan)
% {FC1, FC2, Fz, F3, F4,  AF3, AF4};

% for a set of sensors around FCz, 64 chan
% {FC1, FC2, FCz, C1, C2, CZ};
% {C3, C4, FC3, FC4, Fz1:4 }
% keyChans = [4, 5, 10, 11, 12, 13, 39, 40, 45, 46, 47, 48, 49, 50];

keyHz = [4 7]; % theta is 4-7Hz.
% keyHz = [10 14]; % alpha is 10-14Hz.

timeWindow = [0 8700]; % N1 peak is ~110ms (but wavelet window is ~500ms)

% which condition do we want to process?
% currently written to loop through all conditions.
% allConditions = {'B1(155)' , 'B2(151)' , 'B3(153)' , ...
%     'B4(55)' , 'B5(51)' , 'B6(54)', 'B7(115)', ...
%     'B8(111)' , 'B9(113)' , 'B10(15)', 'B11(11)', ...
%     'B12(13)'};

% just bins 1-3. 
allConditions = {'B1(' , 'B2(' , 'B3(', ...
    'B4(', 'B5(' , 'B6(', ...
    'B7(', 'B8(', 'B9('};
% fixed values across all figures.
colour_max = 7; % or 3
colour_min = -7; % or -1
colour_max_itc = 0.5;
colour_min_itc = 0;

% global colour scheme
colScheme = 'jet';

% what's the relevant config file called?
ConfigFileName = 'Config_CaitlinCharlotte_v2023_onlyCompleteFiles';
rawDataFolder = 'TF_output_dBwholeWindow';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Current_File_Path = pwd;
addpath('Functions');
ConfigFilePath = [Current_File_Path filesep 'SupportingDocs' filesep ConfigFileName '.xlsx'];
Options = detectImportOptions(ConfigFilePath);

for k = 1:numel(Options.VariableTypes)
    Options.VariableTypes{k} = 'char';
end
DataConfig = table2struct(readtable(ConfigFilePath, Options));
DataConfig = adjustConfigData(DataConfig);

%% identify channels and locations given incomplete chan lists%%%%%%%%%%%%%
NoOfChans = DataConfig.TotalChannels{1};

load('chanlocs.mat');

% get channel numbers
% that is, convert channel indices from channel names.
% have done this assuming that all bins have the same channel location
% structure, but robust to missing values for some bins.
if strcmp(keyChans{1}, 'all')
    keyChanIdx = DataConfig.firstScalp:DataConfig.lastScalp;
else
    for ThisChan = 1:length(keyChans)
        if ~isempty(find(strcmp({chanlocs.labels}, keyChans{ThisChan})==1))
            keyChanIdx(ThisChan) = find(strcmp({chanlocs.labels}, keyChans{ThisChan})==1);
        end
    end
end

%% open up eeglab and initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% and open eeglab to access the EEGlab functions
eeglab;
% shorten variable name
SUB = DataConfig.SUB;

%% goal #1, plot  theta spatial topo during ERP peak (100-200ms) for patients
% initialize the topoplot output variables.

spatial_patient_ersp = NaN(length(SUB),size(allConditions, 1), NoOfChans);
spatial_patient_itc = NaN(length(SUB),size(allConditions, 1), NoOfChans);

for thisPID = 1:length(SUB)
    % loop through each needed file, load it, grab needed data, push it into
    % some holder variable.
    filename = [rawDataFolder filesep SUB{thisPID} '_TFdata.mat'];
    load(filename); % create variable tf_data
    tf_data = outData;
    for thisCND = 1:length(allConditions)
        % want 64 chans, average over timeWindow, average over keyHz
        time_idx = (tf_data.cond(thisCND).times >= timeWindow(1)) & ...
            (tf_data.cond(thisCND).times <= timeWindow(2)) ;
        freq_idx = (tf_data.cond(thisCND).freqs >= keyHz(1)) & ...
            (tf_data.cond(thisCND).freqs <= keyHz(2)) ;
        disp(['processing_cnd_' allConditions{thisCND} '_in_PID_' SUB{thisPID} ])
        for thisChan = 1:NoOfChans
            % find the channel in the input file that matches the desired
            % channel label, and make it conform to "chanlocs" variable.

            match = strcmp({tf_data.cond(thisCND).chan.lbl}, chanlocs(thisChan).labels);
            if any(match)
                % that channel appears in the data set.
                [row, col] = find(match);
                rawChan_idx = max(row, col);
                % found the channel, now lets pull those data.
                spatial_patient_ersp(thisPID, thisCND, thisChan) =  ...
                    mean(tf_data.cond(thisCND).chan(rawChan_idx).ersp(freq_idx,time_idx), 'all', 'omitnan');
                spatial_patient_itc(thisPID, thisCND, thisChan) =  ...
                    mean(tf_data.cond(thisCND).chan(rawChan_idx).itc(freq_idx,time_idx), 'all', 'omitnan');

            else
                % do nothing. It's already initialized as a NaN.
            end

        end % of channel by channel loop
    end % of condition by condition loop

end % of PID loop

% create output folders if needed 
if exist('TF_images', 'dir') == 7
else
    mkdir 'TF_images'
end

save('TF_images/spatial_patient_ersp.mat', 'spatial_patient_ersp');
save('TF_images/spatial_patient_itc.mat', 'spatial_patient_itc');

% plot each condition
% colour_max = max(spatial_patient_ersp(:));
% colour_min = min(spatial_patient_ersp(:));
% colour_max_itc = max(abs(spatial_patient_itc(:)));
% colour_min_itc = min(abs(spatial_patient_itc(:)));

for thisCND = 1:length(allConditions)
    % draw each figure for ersp
    figure;
    % just some wrangling with dimenions to make sure it's the
    dataTopoplot = squeeze(mean(spatial_patient_ersp(:,thisCND,1:NoOfChans),1, 'omitnan'));
    % a quick check to see if the data variable is empty.
    if sum(isnan(dataTopoplot)) == length(dataTopoplot)
        % do nothing. Skip itc/ersp spatial plotting. variable is empty
    else

        topoplot(dataTopoplot(1:DataConfig.TotalChannels{1}), chanlocs(1:DataConfig.TotalChannels{1}));
        % sort out the colorBar
        colorbar;
        if ~isempty(colour_min) || ~isempty(colour_max)
            caxis([colour_min, colour_max]);
        end
        colormap(colScheme);
        % name and save each figure.
        f = gcf;
        f.Units = 'inches';
        f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
        fig_filename = ['TF_images' filesep 'Patient_ersp_' allConditions{thisCND} '_TargetHzTopoplot.png'];
        disp(['Saving topoimage image ' fig_filename]);
        exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
        close(gcf);

        % and output a spreadsheet of these values too.
        spreadsheetData = squeeze(spatial_patient_ersp(:,thisCND,:));
        excelFilename = ['TF_images' filesep 'spatial_patient_ersp.xlsx'];
        writematrix(spreadsheetData, excelFilename, 'Sheet', allConditions{thisCND});
        writecell(SUB', excelFilename, 'Sheet', 'SUBS');

        % draw each figure for itc
        figure;
        dataTopoplot = squeeze(mean(abs(spatial_patient_itc(:,thisCND,1:NoOfChans)),1, 'omitnan'));
        topoplot(dataTopoplot(1:DataConfig.TotalChannels{1}), chanlocs(1:DataConfig.TotalChannels{1}));
        % colorBar
        colorbar;
        colormap(colScheme);
        caxis([colour_min_itc, colour_max_itc]);

        % name and save each figure.
        f = gcf;
        f.Units = 'inches';
        f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
        fig_filename = ['TF_images' filesep 'Patient_itc_' allConditions{thisCND} '_TargetHzTopoplot.png'];
        disp(['Saving topoimage image ' fig_filename]);
        exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
        close(gcf);

        % and output a spreadsheet of these values too.
        spreadsheetData = squeeze(abs(spatial_patient_itc(:,thisCND,:)));
        excelFilename = ['TF_images' filesep 'spatial_patient_itc.xlsx'];
        writematrix(spreadsheetData, excelFilename, 'Sheet', allConditions{thisCND});
        writecell(SUB', excelFilename, 'Sheet', 'SUBS');
    end % end of empty check
end % of thisCND loop


%% goal #2, plot  theta spatial topo during ERP peak (100-200ms) for patients


%% goal 3, grab mean time-freq plot for whole epoch, all freqs for patients
for thisPID = 1:length(SUB)
    % loop through each needed file, load it, grab needed data, push it into
    % some holder variable.
    filename = [rawDataFolder filesep SUB{thisPID} '_TFdata.mat'];
    load(filename); % create variable tf_data
    tf_data = outData;
    for thisCND = 1:length(allConditions)
        % want just the key channel, all times, all freqs.
        disp(['processing_cnd_' allConditions{thisCND} '_in_PID_' SUB{thisPID} ])


        for thisChan_idx = 1:length(keyChanIdx)
            % find the list of desired channels, pick them one by one.
            thisChan = keyChanIdx(thisChan_idx);
            % find the channel in the input file that matches the desired
            % channel label, and make it conform to "chanlocs" variable.
            match = strcmp({tf_data.cond(thisCND).chan.lbl}, chanlocs(thisChan).labels);
            if any(match)
                % that channel appears in the data set.
                [row, col] = find(match);
                rawChan_idx = max(row, col);
                % found the channel, now lets find it and pull the data.
                holder_ersp(thisChan_idx, :, :) = tf_data.cond(thisCND).chan(rawChan_idx).ersp;
                holder_itc(thisChan_idx, :, :) = abs(tf_data.cond(thisCND).chan(rawChan_idx).itc);

            else
                % fill in with an appropriately sized NaN matrix.
                holder_ersp(thisChan_idx, :, :) = NaN(size(tf_data.cond(thisCND).chan(rawChan_idx).ersp));
                holder_itc(thisChan_idx, :, :) = NaN(size(tf_data.cond(thisCND).chan(rawChan_idx).itc));
            end

        end % of channel by channel loop

% holder_ersp (chanbs, freqs, times). [~64, 395, 200]
%         for thisChan_idx = 1:length(keyChanIdx)
%             thisChan = keyChanIdx(thisChan_idx);
%             holder_ersp(thisChan_idx, :, :) = tf_data.cond(thisCND).chan(thisChan).ersp;
%             holder_itc(thisChan_idx, :, :) = abs(tf_data.cond(thisCND).chan(thisChan).itc);
%         end % of channel by channel loop


        % 4d: pid, cnd, time, freq.
        temporal_patient_ersp(thisPID, thisCND, :, :) =  ...
            mean(holder_ersp,1, 'omitnan');
        % 4d: pid, cnd, time, freq.
        temporal_patient_itc(thisPID, thisCND, :, :) =  ...
            mean(holder_itc,1, 'omitnan');
    end % of condition by condition loop
end % of PID loop

save('TF_images/temporal_patient_ersp.mat', 'temporal_patient_ersp');
save('TF_images/temporal_patient_itc.mat', 'temporal_patient_itc');

% plot each condition
% colour_max = max(temporal_patient_ersp(:));
% colour_min = min(temporal_patient_ersp(:));
% colour_max_itc = max(abs(temporal_patient_itc(:)));
% colour_min_itc = min(abs(temporal_patient_itc(:)));
times = tf_data.cond(1).times;
freqs = tf_data.cond(1).freqs;

for thisCND = 1:length(allConditions)
    % draw each figure for ersp
figure;
% 4d: pid, cnd, time, freq.
plotMatrix = smooth2a(squeeze(mean(temporal_patient_ersp(:,thisCND,:,:), 1, 'omitnan')), 5, 5);
contourf(times,freqs, plotMatrix );
% contourf(times,freqs, squeeze(mean(temporal_patient_ersp(:,thisCND,:,:), 1, 'omitnan'))  );
colorbar;
colormap(colScheme);
if ~isempty(colour_min) || ~isempty(colour_max)
    caxis([colour_min, colour_max]);
end
% font properties
ax = gca;
ax.FontSize = 16;
ax.FontWeight = 'bold';
ax.LineWidth = 3;

% name and save each figure.
f = gcf;
f.Units = 'inches';
f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
fig_filename = ['TF_images' filesep 'Patient_ersp_' allConditions{thisCND} '_timefreq.png'];
disp(['Saving timefreq image ' fig_filename]);
exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
close(gcf);

% and output a spreadsheet of these values too. 

% draw each figure for itc
figure;
plotMatrix = smooth2a(squeeze(mean(temporal_patient_itc(:,thisCND,:,:), 1, 'omitnan'))  , 5, 5);
contourf(times,freqs, plotMatrix );
% contourf(times,freqs, squeeze(mean(temporal_patient_itc(:,thisCND,:,:), 1, 'omitnan'))  );
colorbar;
colormap(colScheme);
caxis([colour_min_itc, colour_max_itc]);

% font properties
ax = gca;
ax.FontSize = 16;
ax.FontWeight = 'bold';
ax.LineWidth = 3;

% name and save each figure.
f = gcf;
f.Units = 'inches';
f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
fig_filename = ['TF_images' filesep 'Patient_itc_' allConditions{thisCND} '_timefreq.png'];
disp(['Saving timefreq image ' fig_filename]);
exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
close(gcf);

% % and output a spreadsheet of these values too. 
% spreadsheetData = squeeze(abs(spatial_patient_itc(:,thisCND,:)));
% excelFilename = ['TF_images' filesep 'spatial_patient_itc.xlsx'];
% writematrix(spreadsheetData, excelFilename, 'Sheet', allConditions{thisCND});
end

%% goal 4, grab mean time-freq plot for whole epoch, all freqs for control
