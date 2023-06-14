
% plotting code for Cody's data 010623
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loadData = 0; % 1 = loads constituent data into plotData variable
doPlots = 1; % 1 = draw the plots from plotData
peakFreq = 15; % target frequency.
dataSource = 'RESS_output';

%% enter some basic details into a global structure %%%%%%%%%%%%%%%%%%%%%%%
binTimings = struct;
%
% all these timings are relative to the commencement of the 8700ms
% maintenance period.
binTimings(1).baseline = [-2000, 0];
binTimings(1).measureWindow = [0 78000];
binTimings(1).label = 'B1';
%
binTimings(2).baseline = [-2000, 0];
binTimings(2).measureWindow = [0 78000];
binTimings(2).label = 'B2';
%
binTimings(3).baseline = [-2000, 0];
binTimings(3).measureWindow = [0 78000];
binTimings(3).label = 'B3';
%
binTimings(4).baseline = [-2000, 0];
binTimings(4).measureWindow = [0 78000];
binTimings(4).label = 'B4';
%
binTimings(5).baseline = [-2000, 0];
binTimings(5).measureWindow = [0 78000];
binTimings(5).label = 'B5';
%
binTimings(6).baseline = [-2000, 0];
binTimings(6).measureWindow = [0 78000];
binTimings(6).label = 'B6';
%
binTimings(7).baseline = [-2000, 0];
binTimings(7).measureWindow = [0 78000];
binTimings(7).label = 'B7';
%
binTimings(8).baseline = [-2000, 0];
binTimings(8).measureWindow = [0 78000];
binTimings(8).label = 'B8';
%
binTimings(9).baseline = [-2000, 0];
binTimings(9).measureWindow = [0 78000];
binTimings(9).label = 'B9';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% will need this for permutest function. 
addpath('Functions');

if loadData == 1

    for thisBin = 1:numel(binTimings)
        % load the data needed.
        binLbl = binTimings(thisBin).label;
        baseline = binTimings(thisBin).baseline;
        measureWindow = binTimings(thisBin).measureWindow;

        % out_ts is the name of the loaded varialble (just one)
        % structured PID by samples (2D)
        disp(['Reading data from Hilb_' num2str(peakFreq) 'Hz_Cond' binLbl '.mat']);
        load([dataSource filesep 'Hilb_' num2str(peakFreq) 'Hz_Cond' binLbl '.mat']);
        plotData(thisBin).Hilb = out_ts;
        %
        disp(['Reading FFT data from FFT_' num2str(peakFreq) 'Hz_Cond' binLbl '.mat']);
        load([dataSource filesep 'FFT_' num2str(peakFreq) 'Hz_Cond' binLbl '.mat']);
        plotData(thisBin).FFT = out_psd;

        % load in timing data.
        disp(['Reading times from ' num2str(peakFreq) 'Hz_Cond' binLbl '.xlsx']);
        plotData(thisBin).times = readmatrix([dataSource filesep num2str(peakFreq) 'Hz_Cond' binLbl '.xlsx'],'Sheet','times');
        % load in frequency data.
        disp(['Reading hz from ' num2str(peakFreq) 'Hz_Cond' binLbl '.xlsx']);
        plotData(thisBin).hz = readmatrix([dataSource filesep  num2str(peakFreq) 'Hz_Cond' binLbl '.xlsx'],'Sheet','hz');
        % load in PIDs
        disp(['Reading PIDs from ' num2str(peakFreq) 'Hz_Cond' binLbl '.xlsx']);
        plotData(thisBin).PID = readmatrix([dataSource filesep num2str(peakFreq) 'Hz_Cond' binLbl '.xlsx'],'Sheet','SUBs');
        % baseline
        plotData(thisBin).baseline = binTimings(thisBin).baseline;
        % measureWindow
        plotData(thisBin).measureWindow = binTimings(thisBin).measureWindow;
    end

    disp('Saving plotData.mat summary file');
    save ('plotData.mat','plotData');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doPlots == 1

    % load in plotting data in if it's not already there
    if exist('plotData') == 1
        % then the variable is alreadly loaded.
    else % load it in.
        load('plotData.mat');
    end

    clear figInfo % just in case you're rerunning script. Blank slate.

    % load up some info about the particular figure you want drawn.
    figInfo(1).title = 'Critical Trial';
    figInfo(1).bins = {'NoCue', 'Valid', 'Invalid'};
    figInfo(1).lineColours = {[0, 0, 0], [0,1,0], [0,0,1] };
    figInfo(1).shadeColours = {[0.7, 0.7, 0.7], [0.2, 0.8, 0.2], [0.2, 0.2, 0.8]};
    figInfo(1).lineStyles = {':', '-', '-'};
    figInfo(1).binNos = [1,2,3];
    figInfo(1).yLimits = [];
    figInfo(1).useBaseline = 1; % 1 for include baseline, 0 for not.
    figInfo(1).dependMeans = false;
    %
    figInfo(2).title = 'Control Trial';
    figInfo(2).bins = {'NoCue', 'Valid', 'Invalid'};
    figInfo(2).lineColours = {[0, 0, 0], [0,1,0], [0,0,1] };
    figInfo(2).shadeColours =  {[0.7, 0.7, 0.7], [0.2, 0.8, 0.2], [0.2, 0.2, 0.8]};
    figInfo(2).lineStyles = {':', '-', '-'};
    figInfo(2).binNos = [4,5,6];
    figInfo(2).yLimits = [];
    figInfo(2).useBaseline = 1; % 1 for include baseline, 0 for not.
    figInfo(2).dependMeans = false;
    %
    figInfo(3).title = 'Training Trials';
    figInfo(3).bins = {'NoCue', 'Valid', 'Invalid'};
    figInfo(3).lineColours =  {[0, 0, 0], [0,1,0], [0,0,1] };
    figInfo(3).shadeColours =  {[0.7, 0.7, 0.7], [0.2, 0.8, 0.2], [0.2, 0.2, 0.8]};
    figInfo(3).lineStyles = {':', '-', '-'};
    figInfo(3).binNos = [7,8,9];
    figInfo(3).yLimits = [];
    figInfo(3).useBaseline = 1; % 1 for include baseline, 0 for not.
    figInfo(3).dependMeans = false;
    %
    figInfo(4).title = 'ValidityTraining_PermTests';
    figInfo(4).bins = {'Valid', 'Invalid'};
    figInfo(4).lineColours =  {[1, 0, 0], [0,1,0]};
    figInfo(4).shadeColours =  {[0.8, 0.2, 0.2], [0.2, 0.8, 0.2]};
    figInfo(4).lineStyles = {'-', '-'};
    figInfo(4).binNos = [8,9];
    figInfo(4).yLimits = [];
    figInfo(4).useBaseline = 1; % 1 for include baseline, 0 for not.
    figInfo(4).dependMeans = false;
    %
    figInfo(5).title = 'ValidityCriticalTest_PermTests';
    figInfo(5).bins = {'Valid', 'Invalid'};
    figInfo(5).lineColours =  {[1, 0, 0], [0,1,0]};
    figInfo(5).shadeColours =  {[0.8, 0.2, 0.2], [0.2, 0.8, 0.2]};
    figInfo(5).lineStyles = {'-', '-'};
    figInfo(5).binNos = [2,3];
    figInfo(5).yLimits = [];
    figInfo(5).useBaseline = 1; % 1 for include baseline, 0 for not.
    figInfo(5).dependMeans = false;
    %
    figInfo(6).title = 'ValidityControlTest_PermTests';
    figInfo(6).bins = {'Valid', 'Invalid'};
    figInfo(6).lineColours =  {[1, 0, 0], [0,1,0]};
    figInfo(6).shadeColours =  {[0.8, 0.2, 0.2], [0.2, 0.8, 0.2]};
    figInfo(6).lineStyles = {'-', '-'};
    figInfo(6).binNos = [5,6];
    figInfo(6).yLimits = [];
    figInfo(6).useBaseline = 1; % 1 for include baseline, 0 for not.
    figInfo(6).dependMeans = false;

    % and if you want a second figure, put that info here.

    NoOfFigs = numel(figInfo);

    for thisFig = 1:NoOfFigs

        % load in figure-specific info.
        figTitle = figInfo(thisFig).title;
        bins = figInfo(thisFig).bins;
        lineColours = figInfo(thisFig).lineColours;
        shadeColours = figInfo(thisFig).shadeColours;
        lineStyles = figInfo(thisFig).lineStyles;
        binNos = figInfo(thisFig).binNos;
        yLimits = figInfo(thisFig).yLimits;
        useBaseline = figInfo(thisFig).useBaseline;
        dependMeans = figInfo(thisFig).dependMeans;

        figure;
        hold on;
        for k = 1:min([length(bins), length(lineColours), length(shadeColours), length(lineStyles), length(binNos)] )
            %
            thisBin = binNos(k);
            thisBinlbl = bins{k};
            %
            tidx = dsearchn(plotData(thisBin).times, plotData(thisBin).baseline');
            timesToPlot = plotData(thisBin).times(tidx(1):end);
            meanToPlot = mean(plotData(thisBin).Hilb(:,tidx(1):end),1, 'omitnan');
            withinSDdata = plotData(thisBin).Hilb(:,tidx(1):end) - mean(plotData(thisBin).Hilb(:,tidx(1):end),2, 'omitnan');
            %
            SEMToPlot = std(plotData(thisBin).Hilb(:,tidx(1):end),0, 1, 'omitnan')./sqrt(size(plotData(thisBin).Hilb,1));
            % SEMToPlot = std(withinSDdata,0, 1, 'omitnan')./sqrt(size(plotData(thisBin).Hilb,1));
            plus_sem = meanToPlot + SEMToPlot;
            minus_sem = meanToPlot - SEMToPlot;
            %
            timesToPlot = timesToPlot';
            x2 = [timesToPlot, fliplr(timesToPlot)];
            inBetween = [plus_sem, fliplr(minus_sem)];
            %
            fill(x2, inBetween,shadeColours{k}, 'FaceAlpha',0.5, 'LineStyle', 'none');
            %
            line(timesToPlot,meanToPlot, 'Color', lineColours{k}, 'LineWidth', 3, 'LineStyle', lineStyles{k});
            % legend entries
            % if baselines present
            if useBaseline == 1
                if k == 1
                    legVals = {bins{k}, bins{k}, '.' ,'.' ,'.' ,'.' }; % empty entries for baseline lines.
                else
                    legVals = [legVals, {bins{k}, bins{k}, '.', '.' }];
                end
            else
                if k == 1
                    legVals = {bins{k}, bins{k}}; % empty entries for baseline lines.
                else
                    legVals = [legVals, {bins{k}, bins{k}}];
                end
            end


            if useBaseline == 1
                % add baseline info per binType
                xline(plotData(thisBin).baseline(1),  'LineWidth', 1.5,'LineStyle', ':', 'Color', 'k' );
                xline(plotData(thisBin).baseline(2),  'LineWidth', 1.5,'LineStyle', ':', 'Color', 'k' );
                
                if k == 1 % only draw fixed lines once.
                    xline(plotData(thisBin).measureWindow(1),  'LineWidth', 1.5,'LineStyle', '-' );
                    xline(plotData(thisBin).measureWindow(2),  'LineWidth', 1.5,'LineStyle', '-' );
                end
            end

        end % of bin by bin loop.

        % and can insert a permutation test here.
        % must go sample by sample, so loop across samples, and take p-vals and
        % plot.
        if length(figInfo(thisFig).binNos) == 2
            % can do a permutation test comparison between two lines.
            % get the data for those two lines.
            for n = 1:length(binNos)
                thisBin = binNos(n);
                tidx = dsearchn(plotData(thisBin).times, plotData(thisBin).baseline');
                comparisonTimes{n} = plotData(thisBin).times(tidx(1):end);
                comparisonData{n} = plotData(thisBin).Hilb(:,tidx(1):end);
            end
            % find their period of overlap. Limit data, time to that
            % period.
            [C,ia,ib] = intersect(comparisonTimes{1},comparisonTimes{2});
            comparisonData{1} = comparisonData{1}(:,ia);
            comparisonTimes{1} = comparisonTimes{1}(ia);
            comparisonData{2} = comparisonData{2}(:,ib);
            comparisonTimes{2} = comparisonTimes{2}(ib);

            % permutest won't tolerate NaNs.
            % so any subjects with NaN data must go. 
            findNaNs_1 = any(isnan(comparisonData{1}),2);
            findNaNs_2 = any(isnan(comparisonData{2}),2);
            if dependMeans == true
                allNaNs = or(findNaNs_1,findNaNs_2);
                comparisonData{1}(allNaNs,:) = [];
                comparisonData{2}(allNaNs,:) = [];
            else
                comparisonData{1}(findNaNs_1,:) = [];
                comparisonData{2}(findNaNs_2,:) = [];
            end
            
            
            % how many samples do we have? 
            NoOfSamples = length(comparisonTimes{1});

            % initialize an output p-value time series of correct length.
            sigVals = zeros(1,NoOfSamples);

            % using https://au.mathworks.com/matlabcentral/fileexchange/71737-permutest
            disp('Performing cluster-based permutation test. Will take a minute.');

            [clusters, p_values, t_sums, permutation_distribution ] = ...
                permutest(comparisonData{1}',comparisonData{2}',dependMeans,0.05,10000,true,10);
            
            clusterFile = [dataSource filesep 'Clusters_' figTitle '.mat'];
            save(clusterFile, 'clusters' , 'p_values' , 't_sums', 'permutation_distribution');

            NoOfClusters = length(clusters);
            for thisCluster = 1:NoOfClusters
                if p_values(thisCluster) < 0.05
                    sigVals(1,clusters{thisCluster}) = 1;
                    line(comparisonTimes{1},sigVals, 'Color','r' , 'LineStyle', ':', 'LineWidth', 2);
                    legVals = [legVals, {['Cluster ' num2str(thisCluster)]}];
                end
            end

%% Older, non-clustered, simple permutation test. %%%%%%%%%%%%%%%%%%%%%%%%%
%             % critical command for non-cluster permutation test. 
%             disp('performing permutation test')
%             p_vals = NaN(1,NoOfSamples);
%             eff_sizes = NaN(1,NoOfSamples);
%             permutations = 1000;
%             for thisSample = 1:NoOfSamples
%                 [p, observeddifference, effectsize] = ...
%                     permutationTest(comparisonData{1}(:,thisSample), comparisonData{2}(:,thisSample), permutations);
%                 p_vals(thisSample) = p;
%                 eff_sizes(thisSample) = effectsize;
%                 if mod(thisSample, 100) == 0
%                     disp(['Sample number: ' num2str(thisSample) ' of ' num2str(NoOfSamples) ' done.'] );
%                 end
%             end
% 
%             % threshold the p-values to < .05
%             sigVals = p_vals < 0.05; 
%             % now add a line of the appropriate length,  
%             line(comparisonTimes{1},sigVals, 'Color', 'r', 'LineStyle', '-' );
%             % add a legend entry
%             legVals = [legVals, {'significance'}];

        end % of comparison check.
        

        hold off

        % change some formatting and save
        if isempty(yLimits)
        else
            ylim(yLimits);
        end
        lgd = legend(legVals, 'Location', 'eastoutside');
        lgd.FontSize = 8;
        %
        ylabel('Power(SNR)', 'Fontsize', 16);
        xlabel('Time(ms)', 'Fontsize', 16);
        ax = gca;
        ax.FontSize = 16;
        %
        f = gcf;
        f.Units = 'inches';
        f.OuterPosition = [0.5 0.5 7.5 7.5]; % make the figure 7 inches in size.
        fig_filename = [dataSource filesep figTitle '.png'];
        disp('Saving all bin waveform');
        exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
        close(gcf);
    end % of figure by figure loop.


end % of plotting check

