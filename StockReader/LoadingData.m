function [data] = LoadingData(time,tag,database1)
% -------------------------------------------
% Loading stock data from google 
% 
% John Nilsson, 2015-07-14 
% ------------------------------------------

%% -- Initiating script --------------------

% Checking current directory
CurrDir = pwd; 

% Directory of target file 
FileDir = 'C:\Users\John\Desktop\StockReader\dailyDataSuite';
SaveDir = fullfile('C:\Users\John\Desktop\StockReader\DailyData',database1);

%% -- Assigning varaibles ------------------

% --- Assigning time-factors ---
% time.format = 'yyyy/mm/dd';                                 % Time Format 
% time.start = '2008/06/08';                                  % Time Start
% time.end = datestr(datenum(date),time.format);              % Time End     ´

% --- Setting filename ---
file.datename = strrep(time.end, '/', '');                              % Replacing '/' with '_'
file.name = 'StockData';                                                % Setting file name  
file.namecat = strcat(file.name,'_',file.datename,'_',tag,'.csv');      % Concatenating with date

% --- Assigning stocks ---
stock.list = GetStocks(tag);

%% --- Executing Stock Download ---      

cd(FileDir)     % Setting directory of target files

if strcmp(database1 ,'Yahoo')
    % --- Executing extraction command ---
    data = getYahooDailyData(stock.list(:,1), time.start, time.end, time.format);
elseif strcmp(database1 ,'EuroInvest')
    data = getEuroInvestDailyData(stock.list(:,2));
end   
% %% --- Saving to file ---
% cd(SaveDir)                             % Setting save directory 
% save(file.namecat,'data');              % Saving to file
cd(CurrDir)                             % Setting Current directory







