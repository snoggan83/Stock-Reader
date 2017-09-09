function ProcOptimize
% -------------------------------------------
% function ProcOptimize
% Optimizing parameters in ProcessStrat
%
% John Nilsson, 2015-11-21
% ------------------------------------------

% --- Setting Initial Conditions of Solver ---
x0(1) = 22.9944;            % Setting Low RSI 
x0(2) = 90.8223;            % Setting High RSI 
x0(3) = 0.8142;             % Setting procentage stop-loss

x0 =   [27.5871   88.5443    0.7639];

% --- Setting options for optimization routine ---
options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-4);                                           

% --- Executing newton solver ---
[x,fval] = fminsearch(@ProcessStrat,x0,options);                                               

x
