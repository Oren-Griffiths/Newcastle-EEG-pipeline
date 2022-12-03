%% RESS empirical example
% this script  adapted from that which accompanies the manuscript titled
%    "Rhythmic entrainment source separation:
%     Optimizing analyses of neural responses to rhythmic sensory stimulation"
%   from Mike X Cohen and Rasa Gulbinaite
% mikexcohen@gmail.com or rasa.gulbinaite@gmail.com

clear

%% specify parameters
for thisLoop = 1:4
    
    % which experiment are we going to run?
    ConfigFileName = 'Config_Natalie';
    plotting = 0; % 1 = draw figs per person/freq; 0 = don't draw.
    extractData = 1; % 1 = pull raw data; 0 = use existing files.
    

    
    % the experiment had 15Hz, 17.14Hz, 20Hz, 24Hz.
    %   15Hz and 24Hz are targets/distractors.
    %   17.14 and 20Hz are the unexpected object.
    peakFreqs = [15, 17.14, 20, 24]; % hz
    
    % used for 'best-electrode' analyses
    electrode1 = 'Oz';
    
    % condition number, see below
    % nominate some strings that are in the 'eventtype' field for each relevant
    % epoch.
    switch thisLoop
        case 1
            allConds = {'B1'};
        case 2
            allConds = {'B2'};
        case 3
            allConds = {'B3'};
        case 4
            allConds = {'B4'};
    end
    % allConds = {'B4'}; % B4 and B2 done at 10pm.
    baseline = [-16000 -14000];
    measureWindow = [0 14000];
    
    % parameters for RESS filters:
    peakwidt  = .5; % FWHM at peak frequency
    neighfreq = 1;  % distance of neighboring frequencies away from peak frequency, +/- in Hz
    neighwidt = 1;  % FWHM of the neighboring frequencies
    
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
            % and remove any non-scalp data (not needed with Newcastle preprocessing).
            % EEG = pop_select(EEG, 'channel',[ 1:32]);
            
            % do a quick check to see whether "best" electrode was removed in this
            % data set.
            if isempty(find(strcmp({EEG.chanlocs.labels}, electrode1)==1))
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
                index = index(index>0); % remove any empty searches
                
                % check to see that there are actually epochs available.
                if isempty(index)
                else
                    
                    % finally pull out the needed data.
                    data  = EEG.data(:,:,index);
                    
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
                    
                    ressx = mean(abs( fft(ress_ts1(tidx(1):tidx(2),:),nfft,1)/diff(tidx) ).^2,   2    );
                    
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
                    
                    denom_ress = mean(hilbdat_ress(:,EEG.times < baseline(2)),2);
                    % denominator (epochs in dim 1).
                    
                    for thisEpoch = 1:size(hilbdat_ress,1)
                        % report this in epochs by samples format.
                        snr_hilb_ress(thisEpoch, :) = hilbdat_ress(thisEpoch,EEG.times > baseline(2))./denom_ress(thisEpoch);
                    end
                    % average across epochs
                    snr_hilb_ress = mean(snr_hilb_ress,1);
                    times_hilb = EEG.times(EEG.times > baseline(2))';
                    
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
                        topoplot(map2plot./max(map2plot),EEG.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off','shading','interp');
                        title([ 'RESS for ' num2str(peakFreqs(thisFreq)) ' Hz' ])
                        
                        %subplot(243)
                        figure
                        map2plot = dataX(:,dsearchn(hz',peakFreqs(thisFreq)));
                        topoplot(map2plot./max(map2plot),EEG.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off','emarker2',{find(strcmpi({EEG.chanlocs.labels},electrode1)) 'o' 'w' 4},'shading','interp');
                        title([ 'RESS for ' num2str(peakFreqs(thisFreq)) ' Hz' ])
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
                    temp = char(string(allConds));
                    condLbl = temp(:)';
                    outFilename = [Subject_Path filesep ...
                        SUB{thisSUB} '_' num2str(peakFreqs(thisFreq)) 'Hz_Cond' ...
                        condLbl '.mat'];
                    % and out it goes.
                    save(outFilename,'snr_hilb_ress','times_hilb', 'ressx', 'snrR' , 'hz', ...
                        'ress_normWeights', 'chanlocs');
                end % of skipping out of a no-epoch data set.
            end % of freq by freq loop.
        end % of subject by subject loop.
    end % of decision not to extract data
    
    % now go subject by subject and generate a time-by-SUB power spreadsheet
    % one per condition and per frequency.
    
    % these variables defined above.
    % peakFreqs = [15, 17.14, 20, 24]; % hz
    % allConds = {'B1' , 'B2' , 'B3' , 'B4' , 'B1B2B3B4'};
    % SUB = DataConfig.SUB;
    
    % load one file to initialise some basic values
    Subject_Path = [fileparts(pwd) filesep SUB{1}];
    openFile = [Subject_Path filesep SUB{1} '_' num2str(peakFreqs(1)) 'Hz_Cond' ...
        allConds{1}  '.mat'];
    load(openFile);
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
            writematrix(out_ts, outFilename, 'Sheet' , ...
                ['Hilb_' num2str(peakFreqs(thisFreq)) 'Hz_Cond' allConds{thisCond} ]);
            writematrix(out_psd, outFilename, 'Sheet' , ...
                ['FFT_' num2str(peakFreqs(thisFreq)) 'Hz_Cond' allConds{thisCond} ]);
            % save as matlab data too.
            outFilename = ['RESS_output' filesep 'Hilb_' num2str(peakFreqs(thisFreq)) 'Hz_Cond' allConds{thisCond} '.mat'];
            save(outFilename, 'out_ts');
            outFilename = ['RESS_output' filesep 'FFT_' num2str(peakFreqs(thisFreq)) 'Hz_Cond' allConds{thisCond} '.mat'];
            save(outFilename, 'out_psd');
        end % of thisFreq loop
    end % of thisCond loop
    
end % of extra loop across all conditions individually

