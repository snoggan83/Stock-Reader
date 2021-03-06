function indicate = ProcessStrat(scale,params)
% -------------------------------------------
% function ProcessData
% Processing data collected from EuroInvest
%
% John Nilsson, 2015-11-21
% ------------------------------------------

%% --- Initializing script ---

% clear all
close all
clc

direct.curr = pwd;                                                          % Checking Current Directory

%% --- Setting properties ---

% --- Database Selection ---
file.database = 'EuroInvest';                                               % Select Database

% --- Assigning Tags ---
%  tag = 'LargeCap';                                                          % Setting tag
%  tag = 'MidCap';                                                             % Setting tag

% --- Assigning time-factors ---
time.format = 'yyyy/mm/dd';                                                 % Time Format
time.start = '2013/01/01';                                                  % Time Start
time.end = datestr(datenum(date),time.format);                              % Time End

% --- Setting filename ---
file.datename = strrep(time.end, '/', '');                                  % Replacing '/' with '_'
file.name = 'StockData';                                                    % Setting file name
file.namecat = strcat(file.name,'_',file.datename,'_',tag,'.csv');          % Concatenating with date

% --- Settting search mode ---
mode = 'search';                                                            % Searching for relevant stocks
mode = 'Optimize';                                                          % Plotting


direct.save = fullfile('C:\Users\John\Desktop\StockReader\DailyData',...    % Setting full Directorz of save repository
    file.database);

pwd

cd('C:\Users\John\Desktop\StockReader')
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
ki = 1;
vars = fieldnames(data);       % fieldnames

for i = 1:size(par.list,1)
    % --- Collecting induvidual datasets ---
    Mat.dataset = getfield(data, vars{i});
    % --- Converting from dataset to double ---
    Mat.double = double(Mat.dataset);

    % --- Setting time period ---
    time.elm = find(datenum(time.start,time.format) < Mat.double(:,1));
    Mat.double_new = Mat.double(time.elm,:);

    % --- Assigning headers ---
    Headers = get(Mat.dataset, 'VarNames');

    var = []; % Initializing structure
    % --- Creating structure ---
    for j = 1:length(Headers)
        var = setfield(var, Headers{j}, Mat.double_new(:,j));
    end


    %% --- Calculating SMA ---

    % --- Adjusting for vector length ---

    % Single Moving average bands
    SMA.values = [200 50 20];

    for k = SMA.values

        time.elm_SMA = [(time.elm(1) - (k-1)):(time.elm(1)-1) time.elm']';

        if time.elm_SMA(1) <= 0
            time_elm_SMAover = find(time.elm_SMA > 0);
            time.elm_SMA = time.elm_SMA(time_elm_SMAover);
        end

        var.Close_SMA = Mat.double(time.elm_SMA,5);

        SMA1 = conv(var.Close_SMA,ones(k,1))/length(ones(k,1));      % Temporary SMA value, full convolution
        SMA2 = SMA1(length(ones(k,1)):end-length(ones(k,1))+1);      % Selecting valid convolution values
        SMA = setfield(SMA,strcat('n',num2str(k)),SMA2);             % Assigning structure

        clear time.elm_SMA var.Close_SMA
    end

    %% --- Calulating relative strength index (RSI) --

    nperiods = 14;                                                  % Number of periods used in RSI calculation

    % --- Using MATLABs built in RSI function --
    rsi.par = rsindex(var.Close,14);

    %% --- Searching for RSI threshold ---

    if nargin ~= 0
        % --- If optization is executed ---
        RSI.THlow = params.parameter1(1);               % Low RSI Threshold
        RSI.THhigh = params.parameter2(1);              % High RSI Threshold
        RSI.Pros = params.parameter3(1);                % Pull out
        RSI.Der = params.parameter4(1);                 % Derivative of 200 index
    else
        % --- If normal script is executed ---
        RSI.THlow = 3.1800;                                 % Low RSI Threshold
        RSI.THhigh = 89.9200;                           % High RSI Threshold
        RSI.Pros = 0.3600;                              % Pull out
        RSI.Der = 0.0960;                               % Derivative of 200 index
    end

    % --- Initializng ---
    rsi.bool = 'start';
    rsi.in = 1;
    rsi.in1 = 1;

    deffo = 1;

    % --- Adating to length of n200 ---
    rsi.par200 = rsi.par(end-length(SMA.n50)+1:end);
    var.Date200 = var.Date(end-length(SMA.n50)+1:end);
    var.Close200 = var.Close(end-length(SMA.n50)+1:end);

    % --- Calculating gradient of n200 --
    SMA.n200grad = gradient(SMA.n50);

    for j = 2:length(rsi.par200)
        if rsi.par200(j-1) >= RSI.THlow  && rsi.par200(j) < RSI.THlow  && SMA.n200grad(j) > RSI.Der && strcmp(rsi.bool,'start')            % Locating starting point
            rsi.nval_buy(rsi.in) = rsi.par200(j);
            rsi.tval_buy(rsi.in) = var.Date200(j);
            rsi.close_buy(rsi.in) = var.Close200(j);
            deffo = var.Close200(j);
            rsi.bool = 'end';
            rsi.in = rsi.in + 1;
        elseif rsi.par200(j-1) < RSI.THhigh && rsi.par200(j) >= RSI.THhigh && strcmp(rsi.bool,'end') || var.Close200(j) < deffo*RSI.Pros  && strcmp(rsi.bool,'end')    % Locating end point
            rsi.nval_sell(rsi.in1) = rsi.par200(j);
            rsi.tval_sell(rsi.in1) = var.Date200(j);
            rsi.close_sell(rsi.in1) = var.Close200(j);
            rsi.in1 = rsi.in1 + 1;
            rsi.bool = 'start';
        end
    end

    %     if isfield(rsi, 'tval_sell')
    %         figure(n); n = n+1;
    %         plot(var.Date,rsi.par,'-')
    %         hold on
    %         plot(rsi.tval_buy,rsi.nval_buy,'mo','MarkerFaceColor','m');                         % Buy signal
    %         hold on
    %         plot(rsi.tval_sell,rsi.nval_sell,'co','MarkerFaceColor','c');                       % Sell signal
    %     end

    %     grid on                                                                                 % Turning grid on
    %     dynamicDateTicks                                                                        % Setting dynamic axis indication
    %
    %     labl(1) = xlabel('Time');                                                               % Setting xlabel axis
    %     labl(2) = title(par.list(i,3));                                                         % Setting title
    %
    %     set([gca labl],'fontsize',14,'fontangle','italic')                                      % Formatting points

    % --- Evaluating earnings of concept
    if isfield(rsi,'tval_sell')
        if length(rsi.tval_buy) ~= length(rsi.tval_sell)
            rsi.tval_buy = rsi.tval_buy(1:end-1);                                               % Removing last elementin time vector
            rsi.close_buy = rsi.close_buy(1:end-1);                                             % Removing last elementin value vector
        end
        for ih = 1:length(rsi.tval_buy)

            % --- Number of stocks ---
            earn.nr = 1/rsi.close_buy(ih);

            % --- Earnings ---
            earn.diff(ih) = earn.nr*(rsi.close_sell(ih) - rsi.close_buy(ih));                   % Earnings
        end
        earn.sum = sum(earn.diff);                                                              % Total Earnings
        totl(ki) = length(earn.diff);                                                           % Total number of transactions
        tot(ki) = earn.sum;
        ki = ki + 1;                                                                            % Incrementing
    end

    if nargin ==0
        if isfield(rsi, 'tval_sell')
            figure(n); n = n+1;
            plot(var.Date,var.Close,'-')
            hold on
            plot(rsi.tval_buy,rsi.close_buy,'mo','MarkerFaceColor','m');                        % Buy signal
            hold on
            plot(rsi.tval_sell,rsi.close_sell,'co','MarkerFaceColor','c');                      % Sell signal
            hold on
            plot(var.Date(end-length(SMA.n200)+1:end),SMA.n200,'r-','linewidth',1)

            grid on                                                                             % Turning grid on
            dynamicDateTicks                                                                    % Setting dynamic axis indication

            labl(1) = xlabel('Time');                                                           % Setting xlabel axis
            labl(2) = title(par.list(i,3));                                                     % Setting title

            text(0.1,0.1,num2str(earn.sum),'Units','normalized','fontsize',20)                  % Adding Earnings Textbox

            set([gca labl],'fontsize',14,'fontangle','italic')                                  % Formatting axis
        end
    end




    %     %% --- Preparing plotting ---
    %
    %     % --- Adding figure ---
    %     fig1 = figure(i);
    %     screensize = get(0, 'Screensize');
    %     set(fig1,'Position',screensize)
    %
    %     % --- Adding the SMA200 Band ---
    %     ax(1) = subplot(2,1,1);
    %     plot(var.Date,var.Close,'-')
    %     hold on
    %     plot(var.Date(end-length(SMA.n200)+1:end),SMA.n200,'r-','linewidth',1)          % SMA200
    %
    %     grid on                                                                         % Turning grid on
    %     dynamicDateTicks                                                                % Setting dynamic axis indication
    %     labl(1) = xlabel('Time');                                                       % Setting xlabel axis
    %     labl(2) = title(par.list(i,3));                                                 % Setting title
    %     set([gca labl],'fontsize',14,'fontangle','italic')                              % Formatting axis

    %     % --- Adding the RSI Screen ---
    %     ax(2) = subplot(2,1,2);
    %     plot(var.Date,rsi.par,'-')
    %     hold on
    %     plot(var.Date,30*ones(length(rsi.par),1),'k-','linewidth',2)                    % Adding lower RSI band (Oversold)
    %     hold on
    %     plot(var.Date,70*ones(length(rsi.par),1),'k-','linewidth',2)                    % Adding lower RSI band (Overbought)
    %     hold on
    %     plot(var.Date,10*ones(length(rsi.par),1),'r-','linewidth',2)                    % Adding buy band
    %
    %     dynamicDateTicks                                                                % Setting dynamic axis indication
    %     labl(1) = xlabel('Time');                                                       % Setting xlabel axis
    %     set([gca labl],'fontsize',14,'fontangle','italic')                              % Formatting axis
    %
    %
    %     linkaxes(ax,'x');                                                               % Linking axes
    %
    %     if i == 10
    %         return
    %     end

    % --- Clearing varaibles ---
    clear rsi earn

end

if exist('tot')
    % --- Summurizing total earnings ---
    indicate = (1000 - sum(tot) + sum(totl)/7);
else
    indicate = 1000;
end

sum(totl)


% 
%     parameter1: 66.0100
%     parameter2: 88.5400
%     parameter3: 0






