%% Algorithmic Trading with MATLAB(R): More Signals
% In <AlgoTradingDemo3.html AlgoTradingDemo3.m> we saw how to add two
% signals together to get improved results using evolutionary learning.  In
% this demo we'll use extend the approach to three signals: MA, RSI, and
% Williams %R.
%
% Copyright 2010, The MathWorks, Inc.
% All rights reserved.
%% Load in some data
% Again we'll import Bund data sampled minutely
load bund1min
testPts = floor(0.8*length(data(:,4)));
step = 30; % 30 minute interval
Bund = data(1:step:testPts,2:end);
BundClose = Bund(:,3); 
BundV = data(testPts+1:step:end,2:end);
BundVClose = BundV(:,3);
annualScaling = sqrt(250*60*11/step);
cost = 0.01;
addpath('gaFiles')

%% WPR performance

wp_ix.range = {1:500};
wfun = @(x) wprFun(x,Bund,annualScaling,cost);
tic
[wp_ix.maxSharpe,wp_ix.param,wp_ix.sh] = parameterSweep(wfun,wp_ix.range);
toc
wpr(Bund,wp_ix.param,annualScaling,cost)
figure
plot(wp_ix.sh)
ylabel('Sharpe''s Ratio')

%% RSI performance

% (Let's find the best perfrorming set of parameters.  In the interest of
% time, I'll set the threshold to 55 (found earlier).
rs.range = {1:300,1:300,55}; % replace 55 by this to do the sweep 50:5:100});

rsfun = @(x) rsiFun(x,BundClose,annualScaling,cost);
tic
[~,rs.param] = parameterSweep(rsfun,rs.range);
toc
rsi(BundClose,rs.param(1:2),rs.param(3),annualScaling,cost)

%% Lead-Lag performance 

ll.seq = [1:20 300:1:400];
ll.ts  = 25:50;
ll.range = {ll.seq,ll.seq,ll.ts};
llfun = @(x) leadlagFun(x,BundClose,annualScaling,cost);

tic
[~,ll.param] = parameterSweep(llfun,range);
toc

leadlag(BundClose(1:ll.param(3):end),ll.param(1),ll.param(2),...
        annualScaling,cost)

return


%% Generate trading signals
N = 2; M = 396; thresh = 55; P = 2; Q = 110;
sma = leadlag(Bund(:,end),N,M,annualScaling,cost);
srs = rsi(Bund(:,end),[15*Q Q],thresh,annualScaling,cost);
swr = wpr(Bund,param,annualScaling,cost);

signals = [sma srs swr];
names = {'MA','RSI','WPR'};
%% Trading signals
% Plot the "state" of the market represented by the signals
figure
ax(1) = subplot(2,1,1); plot(Bund(:,end));
ax(2) = subplot(2,1,2); imagesc(signals')
cmap = colormap([1 0 0; 0 0 1; 0 1 0]);
set(gca,'YTick',1:length(names),'YTickLabel',names);
linkaxes(ax,'x');

%% Generate initial population
% Generate initial population for signals
close all
I = size(signals,2);
pop = initializePopulation(I);
imagesc(pop)
xlabel('Bit Position'); ylabel('Individual in Population')
colormap([1 0 0; 0 1 0]); set(gca,'XTick',1:size(pop,2))
%% Fitness Function
% Objective is to find a target bitstring (minimum value)
type fitness
%%
% Objective function definition
obj = @(pop) fitness(pop,signals,Bund(:,end),annualScaling,cost)
%%
% Evalute objective for population
obj(pop)
%% Solve With Genetic Algorithm
% Find best trading rule and maximum Sharpe ratio (min -Sharpe ratio)
options = gaoptimset('Display','iter','PopulationType','bitstring',...
    'PopulationSize',size(pop,1),...
    'InitialPopulation',pop,...
    'CrossoverFcn', @crossover,...
    'MutationFcn', @mutation,...
    'PlotFcns', @plotRules,...
    'Vectorized','on');

[best,minSh] = ga(obj,size(pop,2),[],[],[],[],[],[],[],options)

%% Evaluate Best Performer
s = tradeSignal(best,signals);
s = (s*2-1); % scale to +/-1
r  = [0; s(1:end-1).*diff(Bund(:,end))-abs(diff(s))*cost/2];
sh = annualScaling*sharpe(r,0);

% Plot results
figure
ax(1) = subplot(2,1,1);
plot(Bund(:,end))
title(['Evolutionary Learning Resutls, Sharpe Ratio = ',num2str(sh,3)])
ax(2) = subplot(2,1,2);
plot([s,cumsum(r)])
legend('Position','Cumulative Return')
title(['Final Return = ',num2str(sum(r),3), ...
    ' (',num2str(sum(r)/Bund(1,end)*100,3),'%)'])
linkaxes(ax,'x');

%%
sma = leadlag(BundV(:,end),N,M,annualScaling,cost);
srs = rsi(BundV(:,end),[P Q],thresh,annualScaling,cost);
swr = wpr(BundV,param,annualScaling,cost);
signals = [sma srs swr];

s = tradeSignal(best,signals);
s = (s*2-1); % scale to +/-1
r  = [0; s(1:end-1).*diff(BundV(:,end))-abs(diff(s))*cost/2];
sh = annualScaling*sharpe(r,0);

% Plot results
figure
ax(1) = subplot(2,1,1);
plot(BundV(:,end))
title(['Evolutionary Learning Resutls, Sharpe Ratio = ',num2str(sh,3)])
ax(2) = subplot(2,1,2);
plot([s,cumsum(r)])
legend('Position','Cumulative Return')
title(['Final Return = ',num2str(sum(r),3), ...
    ' (',num2str(sum(r)/BundV(1,end)*100,3),'%)'])
linkaxes(ax,'x');
%%
% This result isn't as good as the pure moving average case, but it's a
% step in the right direction compared to the MA+RSI case.  Another exercise
% to try is to use this method to combine different signals that capture
% market dynamics better (say a bear, bull, and sideways market) and
% calibrate using the moving training/validation window discussed in demo 3.
%
% But alas, we're moving on to the next demo, which discusses how you can
% speed up you MATLAB code, for you performance junkies out there.  On to
% <AlgoTradingDemo5.html AlgoTrading5.m>.