
% pull out the individual fixation traces
bins = {'B2' ,'B3'};
inputFiles = { ['SupportingDocs' filesep 'Gorilla_Trial_TEST.xlsx' ]    , ...
   ['SupportingDocs' filesep 'Gorilla_Trial_DividedAttention.xlsx' ] };
for i = 1:length(bins)
    % small adjustment to code to loop across both test trials
    thisBin = bins{i};
    thisInputFile = inputFiles{i};
    
    % load in the list of fixation onsets at test
    FixList = readmatrix(thisInputFile);
    FixTotal = size(FixList,1);
    % initialise an output variable
    % 2s pre and 4s post (at 256Hz) + a time 0
    raw_output = NaN(FixTotal, 1537);
    bl_output = NaN(FixTotal, 1537);
    target_output = NaN(FixTotal, 1537);
    
    for k = 1:FixTotal
        %
        disp([num2str(k) ' of ' num2str(FixTotal) ' fixations']);
        ThisOnset = FixList(k,3) - 14000;
        ThisSUB = FixList(k,2);
        Subject_Path = [fileparts(pwd) filesep num2str(ThisSUB)];
        Notice = FixList(k,4);
        CB = FixList(k,5);
        
        % load the relevant file.
        if CB == 1 % IB at 17.14Hz
            dataFile = [num2str(ThisSUB) '_17.14Hz_Cond' thisBin '.mat'];
            targetFile = [num2str(ThisSUB) '_24Hz_Cond' thisBin '.mat'];
            % distFile = [num2str(ThisSUB) '_15Hz_Cond' thisBin '.mat'];
        else % CB = 2 and IB at 20Hz
            dataFile = [num2str(ThisSUB) '_20Hz_Cond' thisBin '.mat'];
            targetFile = [num2str(ThisSUB) '_15Hz_Cond' thisBin '.mat'];
            % distFile = [num2str(ThisSUB) '_24Hz_Cond' thisBin '.mat'];
        end
        
        if exist([Subject_Path filesep dataFile]) == 2
            IB = load([Subject_Path filesep dataFile]);
            Target = load([Subject_Path filesep targetFile]);
            
            % find the relevant portion of data.
            [time_error,time0] = min(abs(IB.times_hilb - ThisOnset));
            raw_output(k, :) = IB.snr_hilb_ress(time0 - 512:time0+1024);
            baseline = mean(IB.snr_hilb_ress(1:time0));
            bl_output(k, :) = IB.snr_hilb_ress(time0 - 512:time0+1024) - baseline;
            %
            target_baseline = mean(Target.snr_hilb_ress(time0 - 256:time0));
            target_output(k, :) = Target.snr_hilb_ress(time0 - 512:time0+1024) - target_baseline;
        else % no file, so update nothing. Already a NaN.
            disp(['No file for ' num2str(k) ' fixation']);
        end
        
    end
    
    %% plot
    
    % first make sure there's some place to send it. 
        if exist('FixByFix_output', 'dir') == 7
    else
        mkdir 'FixByFix_output'
    end
    
    % basic stats on normalized fix-evoked-SSVEP responses
    bl_mean = mean(bl_output,1, 'omitnan');
    bl_count =  nnz(~isnan(bl_output(:,1)));
    tcrit = tinv(0.975, bl_count-1);
    bl_sem = std(bl_output,0,1, 'omitnan')/sqrt(bl_count);
    bl_sem_hi = bl_mean + bl_sem;
    bl_sem_low = bl_mean - bl_sem;
    bl_CI_hi = bl_mean + bl_sem*tcrit;
    bl_CI_low = bl_mean - bl_sem*tcrit;
    
    % separate out by noticers.
    noticers = (FixList(:,4) == 1);
    bl_not = bl_output(noticers,:);
    bl_not_mean = mean(bl_not,1, 'omitnan');
    bl_not_count =  nnz(~isnan(bl_not(:,1)));
    bl_not_sem = std(bl_not,0,1, 'omitnan')/sqrt(bl_not_count);
    bl_not_sem_hi = bl_not_mean + bl_not_sem;
    bl_not_sem_low = bl_not_mean - bl_not_sem;
    bl_not_CI_hi = bl_not_mean + tcrit*bl_not_sem;
    bl_not_CI_low = bl_not_mean - tcrit*bl_not_sem;
    
    % and then the non-noticers
    nons = (FixList(:,4) == 0);
    bl_non = bl_output(nons,:);
    bl_non_mean = mean(bl_non,1, 'omitnan');
    bl_non_count =  nnz(~isnan(bl_non(:,1)));
    bl_non_sem = std(bl_non,0,1, 'omitnan')/sqrt(bl_non_count);
    bl_non_sem_hi = bl_non_mean + bl_non_sem;
    bl_non_sem_low = bl_non_mean - bl_non_sem;
    bl_non_CI_hi = bl_non_mean + tcrit*bl_non_sem;
    bl_non_CI_low = bl_non_mean - tcrit*bl_non_sem;
    
    % get the x axis
    times = linspace(-2,6, length(bl_mean));
    % time to draw
    figure; % first figure. Both groups individually.
    
    hold on
    
    % calculate the shaded region for noticers
    x2 = [times, fliplr(times)];
    inBetween = [bl_not_CI_hi, fliplr(bl_not_CI_low)];
    fill(x2, inBetween,[0.2, 0.2, 0.8], 'FaceAlpha',0.5, 'LineStyle', 'none');
    % and now draw the mean line plot
    line(times, bl_not_mean, ...
        'Color', [0.2, 0.2, 0.8], 'LineWidth', 2, 'LineStyle', '-' );
    
    x2 = [times, fliplr(times)];
    inBetween = [bl_non_CI_hi, fliplr(bl_non_CI_low)];
    fill(x2, inBetween,[0.8, 0.2, 0.2], 'FaceAlpha',0.5, 'LineStyle', 'none');
    % and now draw the mean line plot
    line(times, bl_non_mean, ...
        'Color', [0.8, 0.2, 0.2], 'LineWidth', 2, 'LineStyle', '-' );
    
     % add some guiding lines
    yline(0); % draw the baseline.
    xline(0); % note the time zero
    
    hold off

    ax = gca;
    ax.FontSize = 16;
    %
    f = gcf;
    f.Units = 'inches';
    f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
    fig_filename = ['FixByFix_output' filesep 'SeparatedByNoticer'  thisBin '.png'];
    disp('Saving images');
    exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
    
    figure;
    
    hold on
    
    % calculate the shaded region for overall vals
    x2 = [times, fliplr(times)];
    inBetween = [bl_CI_hi, fliplr(bl_CI_low)];
    fill(x2, inBetween,[0.5, 0.5, 0.5], 'FaceAlpha',0.5, 'LineStyle', 'none');
    % and now draw the mean line plot
    line(times, bl_mean, ...
        'Color', 'black', 'LineWidth', 2, 'LineStyle', '-' );
    
    % add some guiding lines
    yline(0); % draw the baseline.
    xline(0); % note the time zero
    
    hold off
    
    ax = gca;
    ax.FontSize = 16;
    %
    f = gcf;
    f.Units = 'inches';
    f.OuterPosition = [0.5 0.5 5.5 5.5]; % make the figure 5 inches in size.
    fig_filename = ['FixByFix_output' filesep 'AllTogether'  thisBin '.png'];
    disp('Saving images');
    exportgraphics(f,fig_filename,'Resolution',300); % set to 300dpi and save.
    
    
    
%     % draw individual noticers
%     for thisLine = 1:size(bl_not,1)
%             line(times, bl_not(thisLine,:), ...
%         'Color', 'red', 'LineWidth', 0.5, 'LineStyle', '-' );
%     end
%     
%     % draw individual non-noticers
%         for thisLine = 1:size(bl_non,1)
%             line(times, bl_non(thisLine,:), ...
%         'Color', 'black', 'LineWidth', 0.5, 'LineStyle', '-' );
%         end
    
    
%     
% 
%     % lets have a look at impact on Target too. 
%     % basic stats on normalized fix-evoked-SSVEP responses
%     target_mean = mean(target_output,1, 'omitnan');
%     target_count =  nnz(~isnan(target_output(:,1)));
%     target_sem = std(target_output,0,1, 'omitnan')/sqrt(target_count);
%     target_sem_hi = target_mean + target_sem;
%     target_sem_low = target_mean - target_sem;
%     target_CI_hi = target_mean + target_sem*tcrit;
%     target_CI_low = target_mean - target_sem*tcrit;
%     
%     % get the x axis
%     times = linspace(-2,4, length(bl_mean));
%     
%     % draw the figure now.
%     figure;
%     % draw everyone
%     line(times,bl_mean, 'Color', 'black', 'LineStyle', '-' );
%     line(times,bl_sem_hi, 'Color', 'black', 'LineStyle', '--' );
%     line(times,bl_sem_low, 'Color', 'black', 'LineStyle', '--' );
%     line(times,bl_CI_hi, 'Color', 'black', 'LineStyle', ':' );
%     line(times,bl_CI_low, 'Color', 'black', 'LineStyle', ':' );
%     % draw noticers
%     line(times,bl_not_mean, 'Color', 'red', 'LineStyle', '-' );
%     line(times,bl_not_sem_hi, 'Color', 'red', 'LineStyle', '--' );
%     line(times,bl_not_sem_low, 'Color', 'red', 'LineStyle', '--' );
%     % draw the non noticers
%     line(times,bl_non_mean, 'Color', 'blue', 'LineStyle', '-' );
%     line(times,bl_non_sem_hi, 'Color', 'blue', 'LineStyle', '--' );
%     line(times,bl_non_sem_low, 'Color', 'blue', 'LineStyle', '--' );
%     % line
%     yline(0); % draw the baseline.
%     xline(0); % note the time zero
%     
%     % and a figure for impact on target Hz
%     figure;
%     line(times,smooth(target_mean), 'Color', 'black', 'LineStyle', '-' );
%     line(times,smooth(target_sem_hi), 'Color', 'black', 'LineStyle', '--' );
%     line(times,smooth(target_sem_low), 'Color', 'black', 'LineStyle', '--' );
%     line(times,smooth(target_CI_hi), 'Color', 'black', 'LineStyle', ':' );
%     line(times,smooth(target_CI_low), 'Color', 'black', 'LineStyle', ':' );
%     % line
%     yline(0); % draw the baseline.
%     xline(0); % note the time zero
%     % xlim([-1 4]);
    
    %% save
    if exist('FixByFix_output', 'dir') == 7
    else
        mkdir 'FixByFix_output'
    end
    
    % raw output
    saveName = ['ByFixation_Output_' thisBin];
    writematrix(raw_output, ['FixByFix_output' filesep saveName '.xlsx'], ...
        'Sheet', 'raw_output');
    % bl_output
    writematrix(bl_output, ['FixByFix_output' filesep saveName '.xlsx'], ...
        'Sheet', 'bl_output');
    % effect on target
    writematrix(target_output, ['FixByFix_output' filesep saveName '.xlsx'], ...
        'Sheet', 'target_output');
    % save noticers
        writematrix(bl_not, ['FixByFix_output' filesep saveName '.xlsx'], ...
        'Sheet', 'noticers');
    % save non-noticers separately
            writematrix(bl_non, ['FixByFix_output' filesep saveName '.xlsx'], ...
        'Sheet', 'non-noticers');
    
    save(['FixByFix_output' filesep saveName '.mat'], 'bl_output', 'raw_output', ...
        'target_output');
    
end % of bin by bin loop.

