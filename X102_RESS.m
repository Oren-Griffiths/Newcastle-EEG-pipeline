%% RESS empirical example
% this script  adapted from that which accompanies the manuscript titled
%    "Rhythmic entrainment source separation:
%     Optimizing analyses of neural responses to rhythmic sensory stimulation"
%   from Mike X Cohen and Rasa Gulbinaite
% mikexcohen@gmail.com or rasa.gulbinaite@gmail.com

clear

%
%% other fixed parameters

% which experiment are we going to run?
ConfigFileName = 'HPRnet3_wTraining_140623';
plotting = 0; % 1 = draw figs per person/freq; 0 = don't draw.
extractData = 1; % 1 = pull raw data; 0 = use existing files.

% main parameters
peakFreqs = [15]; % hz
RESSchans = {'all'}; % can have a cell array of names, or 'all', or [] = all
BinsToTest = [1:9]; % which "bins"/epochs are you RESS/FFT/Hilb-ing?
customBins = 0; % if customBins = 1, you will need to specify measurement 
% and baseline windows. 

% used for 'best-electrode' analyses
electrode1 = 'Oz';

% define key to pull individual bins from EEG epoch structure.
for k = 1:length(BinsToTest)
    allConds{k} = ['B' num2str(BinsToTest(k))];
end

% parameters for RESS filters:
peakwidt  = 0.5; % FWHM at peak frequency
neighfreq = 1;  % distance of neighboring frequencies away from peak frequency, +/- in Hz
neighwidt = 0.5;  % FWHM of the neighboring frequencies

% parameter for FFT baseline corrections
skipbins =  2; % .2 Hz (skips two bins; minimal smear anticipated)
numbins  = 10+skipbins; %  1 Hz. Again, assumes minimal smear.

%% setup the study level configuration details.
Current_File_Path = pwd;
addpath('Functions');
ConfigFilePath = [Current_File_Path filesep 'SupportingDocs' filesep ConfigFileName '.xlsx'];
Options = detectImportOptions(ConfigFilePath);

for k = 1:numel(Options.VariableTypes)
    Options.VariableTypes{k} = 'char';
end
DataConfig = table2struct(readtable(ConfigFilePath, Options));

DataConfig = adjustConfigData(DataConfig);
SUB = DataConfig.SUB;

% initialize an output structure.
globalData  = struct;
globalData.hz = [];
globalData.times_hilb = [];
globalData.Hilbert = {};
globalData.FFT = {};
globalData.Weights = {};
globalData.Chans = {};
  
%% you have the option of custom measurement (FFT) and baseline (Hilb) periods.  
binTimings = struct;
if customBins == 0
    % import from DataConfig file.
    for k = 1:length(BinsToTest)
        thisBin = BinsToTest(k);
        binTimings(thisBin).baseline = [DataConfig.BaselineMin{1},DataConfig.BaselineMax{1} ];
        binTimings(thisBin).measureWindow = [DataConfig.BaselineMax{1}, DataConfig.EpochMax{1}];
    end
else % specify in the format below.
    
    %
    binTimings(1).baseline = [-2000, 0];
    binTimings(1).measureWindow = [0 78000];
    %
    binTimings(2).baseline = [-2000, 0];
    binTimings(2).measureWindow = [0 78000];
    %
    binTimings(3).baseline = [-2000, 0];
    binTimings(3).measureWindow = [0 78000];
    %
    binTimings(4).baseline = [-2000, 0];
    binTimings(4).measureWindow = [0 78000];
    %
    binTimings(5).baseline = [-2000, 0];
    binTimings(5).measureWindow = [0 78000];
    %
    binTimings(6).baseline = [-2000, 0];
    binTimings(6).measureWindow = [0 78000];
    %
    binTimings(7).baseline = [-2000, 0];
    binTimings(7).measureWindow = [0 78000];
    %
    binTimings(8).baseline = [-2000, 0];
    binTimings(8).measureWindow = [0 78000];
    %
    binTimings(9).baseline = [-2000, 0];
    binTimings(9).measureWindow = [0 78000];
end

%% start the analysis. 
for counter = 1:length(BinsToTest)
    % allows us to test subsets of the bins.
    thisCond = BinsToTest(counter);
    
    % baselining details
    baseline = binTimings(thisCond).baseline;
    measureWindow = binTimings(thisCond).measureWindow;
    
    if extractData == 1
        %% start the subject by subject loop (thisSUB).
        for thisSUB = 1:length(SUB)
            
            %% load in the data
            
            % contains EEG data (in eeglab format) and leadfield information, which is
            % used to create the anatomical estimations.
            eeglab;
            Subject_Path = [fileparts(pwd) filesep SUB{thisSUB}];
            FileToOpen = [SUB{thisSUB} '_ds_addChans_cleanline_asr_lp_refs_event_weighted_epoch_bl_ar.set'];
            
            disp(FileToOpen);
            
            EEG = pop_loadset( 'filename', FileToOpen, 'filepath', Subject_Path);
            % and remove any non-scalp data (not usually needed with Newcastle preprocessing).
            if EEG.nbchan > DataConfig.lastScalp
                EEG = pop_select(EEG, 'channel',[DataConfig.firstScalp:DataConfig.lastScalp]);
            end
            
            % do a quick check to see whether "best" electrode was removed in this
            % data set.
            if isempty(find(strcmp({EEG.chanlocs.labels}, electrode1)==1, 1))
                OzMissing = 1;
            else
                OzMissing = 0;
            end
            
            %% main loops [thisFreq and thisCond].
            % then start the frequency-by-frequency and condition-by-condition loops.
            for thisFreq = 1:length(peakFreqs)
                
                % find which epochs contain which markers
                epochvect = cell(1,length(EEG.trials));
                for i=1:EEG.trials
                    [~,t0] = min(abs(cell2mat(EEG.epoch(i).eventlatency)));
                    % epochvect(i) = EEG.epoch(i).eventtype{t0};
                    epochvect{i} = EEG.epoch(i).eventtype{t0};
                end
                
                %% Start RESS
                
                % FFT parameters
                nfft = ceil( EEG.srate/.1 ); % .1 Hz resolution
                tidx = dsearchn(EEG.times',measureWindow'); % we use data from .5-10 seconds, to avoid using the stimulus transient
                
                % find the relevant epochs.
                if isempty(find(contains(epochvect,allConds{counter})))
                    index = [];
                else
                    index = find(contains(epochvect,allConds{counter}));
                end
                
                % extra code to help manage when you want to aggregate
                % across multiple bins
                %                 for thisCond = 1:length(allConds)
                %                     if isempty(find(contains(epochvect,allConds{thisCond})))
                %                         index_allConds{thisCond} = 0;
                %                     else
                %                         index_allConds{thisCond} = find(contains(epochvect,allConds{thisCond}));
                %                     end
                %                 end
                %                 % combine all the found epochs
                %                 index = index_allConds{1};
                %                 if length(index_allConds) > 1
                %                     for thisEpoch = 2:length(index_allConds)
                %                         index = union(index, index_allConds{thisEpoch});
                %                     end
                %                 end
                %                 index = index(index>0); % remove any empty searches
                
                % check to see that there are actually epochs available.
                if isempty(index)
                    
                    
                    % empty all the reportable vectors.
                    snr_hilb_ress = []; % mean time series
                    snrR = []; % mean FFT
                    ress_normWeights = [];
                    chanlocs = [];
                    hz = [];
                    times_hilb = [];
                    
                else
                    
                                        % finally pull out the needed data.
                    % but filter in only the channels you need.
                    if isempty(RESSchans)
                        % default is all channels included
                        data  = EEG.data(:,:,index);
                    elseif contains(RESSchans, 'all')
                        % you can specifically include all channels
                        data  = EEG.data(:,:,index);
                    else % we need to find the intersection of the two cell arrays.
                        temp = struct2cell(EEG.chanlocs);
                        listOfChans = squeeze(temp(1,1,:));
                        [C,ia,ib] = intersect(listOfChans,RESSchans); 
                        EEG = pop_select( EEG, 'channel',listOfChans(ia));
                        data = EEG.data(:,:,index);
                    end
                    
                    dataX = mean(abs( fft(data(:,tidx(1):tidx(2),:),nfft,2)/diff(tidx) ).^2,3);
                    hz    = linspace(0,EEG.srate,nfft);
                    
                    % compute covariance matrix at peak frequency
                    fdatAt = filterFGx(data,EEG.srate,peakFreqs(thisFreq),peakwidt);
                    fdatAt = reshape( fdatAt(:,tidx(1):tidx(2),:), EEG.nbchan,[] );
                    fdatAt = bsxfun(@minus,fdatAt,mean(fdatAt,2));
                    covAt  = (fdatAt*fdatAt')/diff(tidx);
                    
                    % compute covariance matrix for lower neighbor
                    fdatLo = filterFGx(data,EEG.srate,peakFreqs(thisFreq)+neighfreq,neighwidt);
                    fdatLo = reshape( fdatLo(:,tidx(1):tidx(2),:), EEG.nbchan,[] );
                    fdatLo = bsxfun(@minus,fdatLo,mean(fdatLo,2));
                    covLo  = (fdatLo*fdatLo')/diff(tidx);
                    
                    % compute covariance matrix for upper neighbor
                    fdatHi = filterFGx(data,EEG.srate,peakFreqs(thisFreq)-neighfreq,neighwidt);
                    fdatHi = reshape( fdatHi(:,tidx(1):tidx(2),:), EEG.nbchan,[] );
                    fdatHi = bsxfun(@minus,fdatHi,mean(fdatHi,2));
                    covHi  = (fdatHi*fdatHi')/diff(tidx);
                    
                    % perform generalized eigendecomposition. This is the meat & potatos of RESS
                    [evecs,evals] = eig(covAt,(covHi+covLo)/2);
                    [~,comp2plot] = max(diag(evals)); % find maximum component
                    evecs = bsxfun(@rdivide,evecs,sqrt(sum(evecs.^2,1))); % normalize vectors (not really necessary, but OK)
                    
                    % extract components and force sign
                    % maps = inv(evecs'); % get maps (this is fine for full-rank matrices)
                    maps = covAt * evecs / (evecs' * covAt * evecs); % this works either way
                    [~,idx] = max(abs(maps(:,comp2plot))); % find biggest component
                    maps = maps * sign(maps(idx,comp2plot)); % force to positive sign
                    % maps = abs(maps);
                    
                    % reconstruct RESS component time series
                    ress_ts1 = zeros(EEG.pnts,size(data,3));
                    for ti=1:size(data,3)
                        ress_ts1(:,ti) = evecs(:,comp2plot)'*squeeze(data(:,:,ti));
                    end
                    
                    %% compute SNR spectrum
                    
                    ressx = mean(abs( fft(ress_ts1(tidx(1):tidx(2),:),nfft,1)/diff(tidx) ).^2,   2 , 'omitnan'   );
                    
                    [snrR,snrE] = deal(zeros(size(hz)));
                    % skipbins =  5; % ignore for baseline correction
                    % numbins  = 20+skipbins; % use for baseline correction.
                    
                    % loop over frequencies and compute SNR
                    for hzi=numbins+1:length(hz)-numbins-1
                        numer = ressx(hzi);
                        denom = mean( ressx([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
                        snrR(hzi) = numer./denom;
                    end
                    
                    % and now let's get the power across time (snr, via Hilbert transform).
                    [filtdat,empVals] = filterFGx(ress_ts1',EEG.srate,peakFreqs(thisFreq),peakwidt,false);
                    hilbdat_ress =  abs(hilbert(filtdat')).^2';
                    
                    % now in epochs by samples format.
                    tidx_base = dsearchn(EEG.times',baseline');
                    denom_ress = mean(hilbdat_ress(:,tidx_base(1):tidx_base(2)),2);
                    % denominator (epochs in dim 1).
                    
                    snr_hilb_ress = [];
                    for thisEpoch = 1:size(hilbdat_ress,1)
                        % report this in epochs by samples format.
                        snr_hilb_ress(thisEpoch, :) = hilbdat_ress(thisEpoch,:)./denom_ress(thisEpoch);
                    end
                    
                    % average across epochs
                    snr_hilb_ress = mean(snr_hilb_ress,1, 'omitnan');
                    times_hilb = EEG.times';
                    
                    % and grab some data about spatial distribution of inferred RESS
                    % component.
                    map2plot = maps(:,comp2plot);
                    ress_normWeights = map2plot./max(map2plot);
                    chanlocs = EEG.chanlocs;
                    
                    %% some plotting...
                    if plotting == 1
                        % figure(1), clf
                        xlim = [peakFreqs(thisFreq)-5 peakFreqs(thisFreq)+5];
                        
                        figure;
                        plot(hz,snrR,'ro-','linew',1,'markersize',5,'markerface','w')
                        hold on
                        plot(hz,snrE,'ko-','linew',1,'markersize',5,'markerface','w')
                        set(gca,'xlim',xlim)
                        axis square
                        xlabel('Frequency (Hz)'), ylabel('SNR')
                        legend({'RESS';electrode1})
                        
                        %
                        figure;
                        map2plot = maps(:,comp2plot);
                        ress_normWeights = map2plot./max(map2plot);
                        if isreal(ress_normWeights)
                            topoplot(ress_normWeights,EEG.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off','shading','interp');
                            title([ 'RESS for ' num2str(peakFreqs(thisFreq)) ' Hz' ])
                        else
                            topoplot(abs(ress_normWeights),EEG.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off','shading','interp');
                            title([ 'COMPLEX eigenvalues, RESS for ' num2str(peakFreqs(thisFreq)) ' Hz' ])
                        end
                        
                        
                        %subplot(243)
                        figure
                        map2plot = dataX(:,dsearchn(hz',peakFreqs(thisFreq)));
                        if isreal(map2plot)
                        else
                            map2plot = abs(map2plot);
                        end
                        topoplot(map2plot./max(map2plot),EEG.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off','emarker2',{find(strcmpi({EEG.chanlocs.labels},electrode1)) 'o' 'w' 4},'shading','interp');
                        title([ 'Electrode power at ' num2str(peakFreqs(thisFreq)) ' Hz' ])
                        
                        % plot the filt_hilb for Ress and Oz
                        figure;
                        line(times_hilb, smooth(snr_hilb_ress, 128), 'Color', 'r');
                        % line(times_hilb, smooth(snr_hilb_Oz, 128), 'Color', 'k');
                        xlabel('Times(s)'), ylabel('SNR')
                        legend({'RESS';electrode1})
                        title([ 'Electrode power across epoch' ])
                        
                        close all;
                    end
                    
                    %% and output the saved data series.
                    % some code specific to "all bins chosen" mode
                    %                     temp = char(string(allConds));
                    %                     condLbl = temp(:)';
                    condLbl = allConds{counter};
                    outFilename = [Subject_Path filesep ...
                        SUB{thisSUB} '_' num2str(peakFreqs(thisFreq)) 'Hz_Cond' ...
                        condLbl '.mat'];
                    % and out it goes.
                    save(outFilename,'snr_hilb_ress','times_hilb', 'ressx', 'snrR' , 'hz', ...
                        'ress_normWeights', 'chanlocs');
                    
                end % of skipping out of a no-epoch data set.
                
                % organized by subject and bin.
                globalData.Hilbert{thisSUB, counter} = snr_hilb_ress; % mean time series
                globalData.FFT{thisSUB, counter} = snrR; % mean FFT
                
                % grab location data per person, if not already grabbed.
                if  length(globalData.Weights) < thisSUB
                    globalData.Weights{thisSUB} = ress_normWeights;
                    globalData.Chans{thisSUB} = chanlocs;
                    if isempty(globalData.Weights{thisSUB})
                        globalData.Weights{thisSUB} = ress_normWeights;
                        globalData.Chans{thisSUB} = chanlocs;
                    end
                end
                
                % get time and hz information, if needed. 
                if isempty(globalData.hz) 
                    globalData.hz = hz;
                    globalData.times_hilb = times_hilb;
                end
                
            end % of freq by freq loop.
        end % of subject by subject loop.
    end % of decision not to extract data
end% of condition by condition loop

% now go subject by subject and generate a time-by-SUB power spreadsheet
% one per condition and per frequency.


% load one file to initialise some basic values
Subject_Path = [fileparts(pwd) filesep SUB{1}];

openFile = [Subject_Path filesep SUB{1} '_' num2str(peakFreqs(1)) 'Hz_Cond' ...
    allConds{1}  '.mat'];

% this should be a while loop...
openFile = [Subject_Path filesep '1_15Hz_CondB7.mat'];
if exist(openFile) > 0
    load(openFile);
else
    disp('Missing Initial file.');
end
%
pnts_out = length(times_hilb);
times_out = times_hilb;
hz_out = hz;
SUB_out = SUB';

% initialize an output array
out_ts= NaN(length(SUB), pnts_out);
out_psd = NaN(length(SUB), length(hz_out));

% initialize an output folder and file
if exist('RESS_output', 'dir') == 7
else
    mkdir 'RESS_output'
end

for thisCond = 1:length(allConds)
    for thisFreq = 1:length(peakFreqs)
        out_ts = NaN(length(SUB),length(times_hilb));
        out_psd = NaN(length(SUB),length(hz));
        
        for thisSUB = 1:length(SUB)
            disp(['loading SUB' SUB{thisSUB}]);
            % load a preprocessed file.
            Subject_Path = [fileparts(pwd) filesep SUB{thisSUB}];
            openFile = [Subject_Path filesep SUB{thisSUB} '_' num2str(peakFreqs(thisFreq)) ...
                'Hz_Cond' allConds{thisCond}  '.mat'];
            if exist(openFile) == 2
                load(openFile);
                % load in the hilb data.
                out_ts(thisSUB,:) = snr_hilb_ress';
                out_psd(thisSUB,:) = snrR;
            end
        end % of thisSUB loop
        
        % now export that as an .xlsx file
        outFilename = ['RESS_output' filesep num2str(peakFreqs(thisFreq)) 'Hz_Cond' allConds{thisCond} '.xlsx'];
        display(['writing file' num2str(peakFreqs(thisFreq)) 'Hz_Cond' allConds{thisCond}] );
        writecell(SUB_out,outFilename, 'Sheet', 'SUBs');
        writematrix(times_out, outFilename, 'Sheet', 'times');
        writematrix(hz_out, outFilename, 'Sheet', 'hz');
        writematrix(out_ts', outFilename, 'Sheet' , ...
            ['Hilb_' num2str(peakFreqs(thisFreq)) 'Hz_Cond' allConds{thisCond} ]);
        writematrix(out_psd, outFilename, 'Sheet' , ...
            ['FFT_' num2str(peakFreqs(thisFreq)) 'Hz_Cond' allConds{thisCond} ]);
        % save as matlab data too.
        outFilename = ['RESS_output' filesep 'Hilb_' num2str(peakFreqs(thisFreq)) 'Hz_Cond' allConds{thisCond} '.mat'];
        save(outFilename, 'out_ts');
        outFilename = ['RESS_output' filesep 'FFT_' num2str(peakFreqs(thisFreq)) 'Hz_Cond' allConds{thisCond} '.mat'];
        save(outFilename, 'out_psd');
    end % of thisFreq loop
end % of condition by condition loop

