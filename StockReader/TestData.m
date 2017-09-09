function TestData
% -------------------------------------------
% function TestData
% Testing data collected from EuroInvest
% Back-Testing Edge Concepts
%
% John Nilsson, 2017-04-08
% ------------------------------------------

%% --- Initializing script ---

clear all
close all
clc

direct.curr = pwd;                                                          % Checking Current Directory
%% --- Setting properties ---

% --- Database Selection ---
file.database = 'EuroInvest';                                               % Select Database

% --- Assigning Tags ---
% tag = 'LargeCap';                                                           % Setting tag
% tag = 'MidCap';                                                           % Setting tag
% tag = 'SmallCap';                                                         % Setting tag
tag = 'Mylist';                                                           % Setting tag
% tag = 'Commodities';                                                      % Setting tag
% tag = 'Index';                                                            % Setting tag

% --- Assigning time-factors ---
time.format = 'yyyy/mm/dd';                                                 % Time Format
time.start = '2010/01/01';                                                  % Time Start
time.end = datestr(datenum(date),time.format);                              % Time End

% --- Setting filename ---
file.datename = strrep(time.end, '/', '');                                  % Replacing '/' with '_'
file.name = 'StockData';                                                    % Setting file name
file.namecat = strcat(file.name,'_',file.datename,'_',tag,'.csv');          % Concatenating with date

% --- Settting search mode ---
% mode = 'search';                                                          % Searching for Relevant stocks
mode = 'plot';                                                            % Plotting

direct.save = fullfile('C:\Users\John\Desktop\StockReader\DailyData',...    % Setting full Directorz of save repository
    file.database);

par.list = GetStocks(tag);                                                  % Assigning list of stocks

%% --- Loading data ---

cd(direct.save)
if exist(file.namecat)
    data = importdata(file.namecat);
    cd(direct.curr)
else
    cd(direct.curr)
    data = LoadingData(time,tag,file.database);
end

%% --- Assigning data ---

n = 1;
for i = 1:size(par.list,1)
    vars = fieldnames(data);       % fieldnames
    
    % --- Collecting induvidual datasets ---
    Mat.dataset = getfield(data, vars{i});
    
    % --- Converting from dataset to double ---
    Mat.double = table2array(Mat.dataset);
    
    % --- Setting time period ---
    time.elm = find(datenum(time.start,time.format) < Mat.double(:,1));
    Mat.double = Mat.double(time.elm,:);
    
    % --- Assigning headers ---
    Headers = Mat.dataset.Properties.VariableNames;
    
    var = []; % Initializing structure
    % --- Creating structure ---
    for j = 1:length(Headers)
        var = setfield(var, Headers{j}, Mat.double(:,j));
    end
    
    
    %% --- Indicators ---
    
    % --- SMA (Trend) ---
    SMA.values = [200 50 20];
    for k = SMA.values
        SMA_temp = tsmovavg(var.Close','s',k);
        SMA = setfield(SMA,strcat('n',num2str(k)),SMA_temp);
    end
    
    % --- MACD (Trend) ---
    [var.macdvec,var.nineperma] = macd(var.Close);
    
    % --- RSI (Momentum) --
    nperiods = 14;
    % --- Using MATLABs built in function --
    var.rsi = rsindex(var.Close,14);
    
    % ADX (Trend) (Avergage Directional Index)
    adx = indicators([var.High,var.Low,var.Close],'adx',14);
    
    % ATR (Volatility) (Avergage Tru Range)
    atr = indicators([var.High,var.Low,var.Close],'atr',14);
    
    DMI_delta = adx(:,1)-adx(:,2);
    idx1 = find(DMI_delta > 0);
    
    
    %% --- Back Testing ---
    
    rsi_low = 25;
    rsi_high = 80;
    
    [s, r, ri, SMA, sh] = rsi(var.Close,14,[rsi_low rsi_high]);
    
    %% Plotting
    
    obv = onbalvol(var.Close,var.Volume);
    
    figure(n), n = n + 1;
    candle(var.High, var.Low, var.Close, var.Open,'k',var.Date,'yy-mm')
    ch = get(gca,'children');
    set(ch(1),'FaceColor','r')
    set(ch(2),'FaceColor','g')
    set(ch(1:2),'EdgeColor','none')
    
    grid on
    labl = title(par.list(i,3));       
    
    
    % figure(2)
    % plot([s,ri]), grid on
    % hold on
    % plot(ones(size(s))*[rsi_low rsi_high],'r.','linewidth',2)
    % legend('Position','Cumulative Return')
    %
    % figure(3)
    % plot(var.Date,[s,cumsum(r)]), grid on
    % legend('Position','Cumulative Return',2)
    % title(['Final Return = ',num2str(sum(r),3),' (',num2str(sum(r)/var.Close(1)*100,3),'%)'])
    
    
    
    figure(n) , n = n + 1;
    ax1 = subplot(4,1,1);
    plot(var.Date,var.Close), grid on
    hold on
    plot (var.Date,SMA,'r-')
    % legend('Position','Cumulative Return',3)
    title(['Final Return = ',num2str(sum(r),3),' (',num2str(sum(r)/var.Close(1)*100,3),'%)'])
    % hold on
    % plot(var.Date(s),var.Close(s),'r.')
    ax2  = subplot(4,1,2);
    
    bar(var.Date,DMI_delta,'facecolor','r');
    hold on
    bar(var.Date(idx1),DMI_delta(idx1),'facecolor','g','edgecolor','g');
    
    % plot(var.Date,adx(:,1),'g-',var.Date,adx(:,2),'r-')
    ax3 = subplot(4,1,3);
    plot(var.Date,ri)
    hold on
    plot(var.Date,ones(size(s))*[rsi_low rsi_high],'r-','linewidth',2);
    ax4 = subplot(4,1,4);
    plot(var.Date,atr)
    linkaxes([ax4,ax3,ax2,ax1],'x');
    
    labl = title(par.list(i,3));                                                     % Setting title
    
    
    % figure(5)
    % % plot(var.Date,obv), grid on
    % % hold on
    % plot(gradient(var.Close),var.Volume,'g*'), grid on
    
end
    
    
    
    
    
    
    
    
    
