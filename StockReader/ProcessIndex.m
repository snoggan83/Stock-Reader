function ProcessIndex
% -------------------------------------------
% function ProcessIndex
% Processing Index data collected from EuroInvest
%
% John Nilsson, 2015-12-30
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

tag = 'Index';                                                              % Setting tag
%tag = 'Commodities';                                                        % Setting tag

% --- Assigning time-factors ---
time.format = 'yyyy/mm/dd';                                                 % Time Format
time.start = '2011/01/01';                                                  % Time Start
time.end = datestr(datenum(date),time.format);                              % Time End

% --- Setting filename ---
file.datename = strrep(time.end, '/', '');                                  % Replacing '/' with '_'
file.name = 'IndexData';                                                    % Setting file name
file.namecat = strcat(file.name,'_',file.datename,'_',tag,'.csv');          % Concatenating with date

% --- Settting search mode ---
mode = 'search';                                                            % Searching for relevant stocks
% mode = 'plot';                                                            % Plotting


direct.save = fullfile('C:\Users\John\Desktop\StockReader\DailyData',...    % Setting full Directorz of save repository
    file.database);

par.list = GetStocks(tag);                                                  % Assigning list of Index

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
vars = fieldnames(data);       % fieldnames

for i = 1:size(par.list,1)
    % --- Collecting induvidual datasets ---
    Mat.dataset = getfield(data, vars{i});
    % --- Converting from dataset to double ---
    Mat.double = double(Mat.dataset);

    % --- Setting time period ---
    time.elm = find(datenum(time.start,time.format) < Mat.double(:,1));
    Mat.double = Mat.double(time.elm,:);

    % --- Assigning headers ---
    Headers = get(Mat.dataset, 'VarNames');

    var = []; % Initializing structure
    % --- Creating structure ---
    for j = 1:length(Headers)
        var = setfield(var, Headers{j}, Mat.double(:,j));
    end

    %% --- Calulating SMA ---

    % ---------------------------------------------------------------------
    % Calculatng Single Moving averages

    %"Golden Cross": A term for a 50-day moving average crosses the 200-day average in an upward motion.
    %Indicates: Trend reversal pointing to a positive trend

    %"Death Cross": A term for a 50-day moving average crosses the 200-day average in a downward motion.
    %Indicates: Trend reversal pointing to a negative trend.
    % ---------------------------------------------------------------------

    % Single Moving average bands
    SMA.values = [200 50 20];

    for k = SMA.values
        SMA1 = conv(var.Close,ones(k,1))/length(ones(k,1));          % Temporary SMA value, full convolution
        SMA2 = SMA1(length(ones(k,1)):end-length(ones(k,1))+1);      % Selecting valid convolution values
        SMA = setfield(SMA,strcat('n',num2str(k)),SMA2);             % Assigning structure
    end

    %% --- Calulating EMA ---

    % Exponential Moving average bands
    EMA.values = [26 12];

    try
        % --- Calculating EMA --
        EMAvec = EMAcalc(var.Close,EMA.values);
        notwork = 'false';
    catch
        notwork = 'true';
        continue
    end

    % --- Setting structure ---
    for EMAind = 1:length(EMAvec)
        EMA = setfield(EMA,strcat('n',num2str(EMA.values(EMAind))),EMAvec{EMAind});
    end

    EMA.grad = gradient(EMA.n12) ;                         % Calculating gradient
    EMA.grad1 = EMA.grad(end-length(EMA.n26)+1:end);       % Modifying gradient-vector
    EMA.t26 = var.Date(end-length(EMA.n26)+1:end);         % Modifying time-vector

    l_up = 1;
    l_down = 1;
    bool = [];
    % --- Locating cross-overs ---
    res = EMA.n26 - EMA.n12(end-length(EMA.n26)+1:end);
    for h = 2:length(res)
        %         if strcmp(bool,'pos') && EMA.grad1(h) < -0.1
        %             EMA.tcross_down(l_down) = EMA.t26(h);
        %             EMA.ncross_down(l_down) = EMA.n26(h);
        %             bool = 'neg';
        %             l_down = l_down+1;
        if res(h) > 0 && res(h-1) < 0
            EMA.tcross_down(l_down) = EMA.t26(h);
            EMA.ncross_down(l_down) = EMA.n26(h);
            bool = 'neg';
            l_down = l_down+1;
        elseif res(h) < 0 && res(h-1) > 0
            EMA.tcross_up(l_up) = EMA.t26(h);
            EMA.ncross_up(l_up) = EMA.n26(h);
            bool = 'pos';
            l_up = l_up+1;
        end
    end

    %% --- Calulating MACD ---

    finobj = [];

    try
        % --- Calculating the EMA9 from the MACD line ----
        EMA.nMACD9 = EMAcalc(-res,9);
    catch
        continue
    end

    EMA.tMACD9 = EMA.t26(end-length(EMA.nMACD9{1})+1:end);

    %% --- Calculating MACD histogram ---

    % --- MACD hist ---
    resMACD = -res(end-length(EMA.tMACD9)+1:end) - EMA.nMACD9{1};

    %% --- Locating turning point ---

    MACD.thold = 0;                                % Threshold limit
    % --- Locating points below threshold ---

    MACD.grad = diff(resMACD);                      % Calculating gradient of Affected points

    ij = 1;
    ij1 = 1;
    for elmloc = 2:length(MACD.grad)-2
        if MACD.grad(elmloc-1) < 0 &&  MACD.grad(elmloc) < 0 && MACD.grad(elmloc+1) > 0  && MACD.grad(elmloc+2) > 0 && resMACD(elmloc) < MACD.thold   % && MACD.grad(elmloc+2) > 0
            elmloc1(ij) =  elmloc;
            for ik = elmloc1(ij)+1:length(MACD.grad)
                if MACD.grad(ik) < 0
                    elmloc2(ij1) = ik;
                    ij1 = ij1 + 1;
                    break
                end
            end

            ij = ij+1;
        end
    end

    %     exist(num2str(elmloc1))
    %     % --- Checking if empty ---
    if ~exist('elmloc1')
        continue
    end
    emptyloc = isempty(elmloc1);

    %% --- Calulating relative strength index (RSI) --

    nperiods = 100;                                              % Number of periods used in RSI calculation

    % --- Using MATLABs built in function --
    rsi.par = rsindex(var.Close,nperiods);

    %% --- Bollinger bands ---

    % Calulaitng bolliger bands
    [boll.mid, boll.uppr, boll.lowr] = bollinger(var.Close);                       % Bollinger Bands with finacial toolbox

    %% --- Keltner Channels ---

    % --- Finding elements of NaN ---
    kelt.High = find(isnan(var.High));
    kelt.Low = find(isnan(var.Low));
    % --- Replacing NaN elements with closing values ---
    var.High(kelt.High) = var.Close(kelt.High);
    var.Low(kelt.Low) = var.Close(kelt.Low);

    % Calculating keltner channels
    vout = indicators([var.High,var.Low,var.Close],'keltner');

    kelt.mid = vout(:,1);
    kelt.uppr = vout(:,2);
    kelt.lowr = vout(:,3);

    %% --- Locating Bollinger squeeze ---

    sqz.in = find(boll.uppr <= kelt.uppr & boll.lowr >= kelt.lowr);                % Locating elemets when both bollinger band are within the keltner channels
    sqz.out = find(boll.uppr >= kelt.uppr & boll.lowr >= kelt.lowr);               % Locating elemets when both bollinger band are outside the keltner channels

    % --- Calculating momentum ---
    mom.val = tsmom(var.Close);                                                    % Calculating momentum

    mom.x = var.Date;                                                              % Calculating x-axis
    mom.y = zeros(length(var.Date),1);                                             % Calculating y-axis

    %% --- Preparing for plotting ---

    indicat = EMA.tMACD9(elmloc1+3);

    %% --- Adding time series ---
    fig1 = figure(n); n = n+1;
    screensize = get(0,'Screensize');
    set(fig1,'Position',screensize)

    if strcmp(tag,'Commodities')
        plot(var.Date,var.Close,'k-.','linewidth',2)                                    % Due to missing information about open, high and low prices
    else
        candle(var.High, var.Low, var.Close, var.Open,'k',var.Date)
        ch = get(gca,'children');
        set(ch(1),'FaceColor','r')
        set(ch(2),'FaceColor','g')
        set(ch(1:2),'EdgeColor','none')
    end

    hold on
    plot(var.Date(end-length(SMA.n200)+1:end),SMA.n200,'r-','linewidth',1)              % SMA200
    hold on
    plot(var.Date(end-length(SMA.n50)+1:end),SMA.n50,'g-','linewidth',1)                % SMA50
    hold on
    plot(var.Date(end-length(SMA.n20)+1:end),SMA.n20,'c-','linewidth',1)                % SMA20
    hold off

    grid on                                                                             % Turning grid on
    dynamicDateTicks                                                                    % Setting dynamic axis indication
    labl(1) = xlabel('Time');                                                           % Setting xlabel axis
    labl(2) = ylabel('Price [SEK]');                                                    % Setting ylabel axis
    labl(3) = title(par.list(i,3));                                                     % Setting title
    set([gca labl],'fontsize',14,'fontangle','italic')                                  % Formatting axes ant titles
    set(gca,'ylim',[min(var.Close) max(var.Close)])                                     % Setting min and max bounds

    %% --- Adding MACD ---
    fig1 = figure(n); n = n+1;
    screensize = get(0,'Screensize');
    set(fig1,'Position',screensize)

    MACD.colelm = find(resMACD > 0);                                                    % Locating values above 0
    MACDhist = bar(EMA.tMACD9,resMACD,'facecolor','r');                                 % Evaluating bar plot
    hold on
    if ~isempty(MACD.colelm)
        MACDhist = bar(EMA.tMACD9(MACD.colelm),resMACD(MACD.colelm)...
            ,'facecolor','g','edgecolor','g');  % Evaluating bar plot
        hold on
    end
    alpha(0.5)                                                                      % Setting transperancy alpha

    % --- Adding Squeeze-indicator ---
    %           plot(mom.x,mom.y,'yo','markerfacecolor','y')
    %           hold on
    plot(mom.x(sqz.in),mom.y(sqz.in),'yo','markerfacecolor','y')                    % Adding squeeze indicator
    %           hold on
    %           plot(mom.x(sqz.out),mom.y(sqz.out),'go','markerfacecolor','g')
    hold off

    grid on                                                                             % Turning grid on
    dynamicDateTicks                                                                    % Setting dynamic axis indication
    labl(1) = xlabel('Time');                                                           % Setting xlabel axis
    labl(2) = title(par.list(i,3));                                                     % Setting title
    set([gca labl],'fontsize',14,'fontangle','italic')                                  % Formatting axes ant titles
    set(gca,'ylim',[min(resMACD) max(resMACD)])                                         % Setting min and max bounds

    %% --- Adding RSI ---
    fig1 = figure(n); n = n+1;
    screensize = get(0,'Screensize');

    set(fig1,'Position',screensize)
    plot(var.Date,rsi.par,'-')
    hold on
    plot(var.Date,30*ones(length(rsi.par),1),'k-','linewidth',2)                    % Adding lower RSI band (Oversold)
    hold on
    plot(var.Date,70*ones(length(rsi.par),1),'k-','linewidth',2)                    % Adding lower RSI band (Overbought)
    hold on
    plot(var.Date,10*ones(length(rsi.par),1),'r-','linewidth',2)                    % Adding buy band

    if strcmp(tag,'Mylist')
        hold on
        plot(var.Date(elmo),rsi.par(elmo),'ro','markerfacecolor','r')                % Plotting starting-point
    end

    grid on                                                                         % Turning grid on
    dynamicDateTicks                                                                % Setting dynamic axis indication
    labl(1) = xlabel('Time');                                                       % Setting xlabel axis
    labl(2) = title(par.list(i,3));                                                 % Setting title
    set([gca labl],'fontsize',14,'fontangle','italic')

    clear EMA MACD1 elmloc1 elmloc2

%     %% --- Identifying candlestick patterns ---
% 
%     % 3 line strike,
% 
%     im = 1;
%     for ij = 4:length(var.Close)
%         if all((var.Close(ij-3:ij-1) - var.Open(ij-3:ij-1)) < 0) && var.Open(ij-1) - var.Open(ij-2) < 0  && var.Open(ij-2) - var.Open(ij-3) < 0 ...
%                 && var.Close(ij-1) - var.Close(ij-2) < 0  && var.Close(ij-2) - var.Close(ij-3) < 0 && (var.Close(ij) - var.Open(ij)) >  (var.Open(ij-3) - var.Close(ij-1)) ...
%                 && var.Open(ij) < var.Close(ij-1) 
%             elm(im) = ij;
%             im = im + 1;
%         end
%     end
% 
%     if exist('elm')
%     
%     fig1 = figure(n); n = n+1;
% 
%     candle(var.High, var.Low, var.Close, var.Open,'k',var.Date)
%     ch = get(gca,'children');
%     set(ch(1),'FaceColor','r')
%     set(ch(2),'FaceColor','g')
%     set(ch(1:2),'EdgeColor','none')
% 
%     hold on 
%     plot(var.Date(elm),var.Close(elm),'o')
%     
%     
%     grid on                                                                             % Turning grid on
%     dynamicDateTicks                                                                    % Setting dynamic axis indication
%     labl(1) = xlabel('Time');                                                           % Setting xlabel axis
%     labl(2) = ylabel('Price [SEK]');                                                    % Setting ylabel axis
%     labl(3) = title(par.list(i,3));                                                     % Setting title
%     set([gca labl],'fontsize',14,'fontangle','italic')                                  % Formatting axes ant titles
%     set(gca,'ylim',[min(var.Close) max(var.Close)])                                     % Setting min and max bounds
%     end
%     return
end











