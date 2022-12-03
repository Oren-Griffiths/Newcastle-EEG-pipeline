
% goes through and checks epoch counts per participant

%% first, setup the config parameters
ConfigFileName = 'Config_Natalie';

% how many bins should there be?
NoOfBins = 4;

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

%% Now loop through and open and review each file.

% intitialize an output variable
epochCounts = NaN(length(SUB), NoOfBins+1);

for thisSUB = 1:length(SUB)
    % creats a variable called "GoodTrials"
    % it is a struct with fields "data" and "ID" (and others)
    Subject_Path = [fileparts(pwd) filesep SUB{thisSUB} ];
    fileName = [Subject_Path filesep SUB{thisSUB} '_ARcorrectedBins.mat'];
    load(fileName);
    disp(fileName);
    
    % no guarantee bins are in intended order, or don't have gaps etc.
    for thisRow = 1:numel(GoodTrials)
        epochCounts( thisSUB  , 1) = str2num(SUB{thisSUB});
        if size(size(GoodTrials(thisRow).data)) < 3
            % if data array is 2D, then just one epoch.
          epochCounts( thisSUB  , GoodTrials(thisRow).ID+1) = 1; 
        else
            % if data array is 3D, then > 1 epochs. 
            epochCounts( thisSUB  , GoodTrials(thisRow).ID+1) = ...
                size(GoodTrials(thisRow).data,3);
        end
    end
end

writematrix(epochCounts, 'epochCounts.xlsx');