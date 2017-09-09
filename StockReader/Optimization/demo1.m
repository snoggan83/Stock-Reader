function demo1
%DEMO1  Demo for usage of DIFFERENTIALEVOLUTION.
%   DEMO1 starts searching the minimum of Rosenbrock's saddle as a demo.
%   Modify this function for your first optimization.
%
%   <a href="differentialevolution.html">differentialevolution.html</a>
%   <a href="http://www.mathworks.com/matlabcentral/fileexchange/18593">File Exchange</a>
%   <a href="https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=KAECWD2H7EJFN">Donate via PayPal</a>
%
%   Markus Buehren
%   Last modified 31.08.2014
%
%   See also DIFFERENTIALEVOLUTION, ROSENBROCKSADDLE.

% set title
optimInfo.title = 'Demo 1 (Rosenbrock''s saddle)';

% --- Checking current directory --- 
curr_dir = pwd;

% --- Setting directory of objective function ---
cd('C:\Users\John\Desktop\StockReader')

% specify objective function
objFctHandle = @ProcessStrat;

cd(curr_dir)

% define parameter names, ranges and quantization:

% 1. column: parameter names
% 2. column: parameter ranges
% 3. column: parameter quantizations
% 4. column: initial values (optional)

paramDefCell = {
	'parameter1', [0 30], 0.01
	'parameter2', [0 100], 0.01
    'parameter3', [0 1], 0.01
    'parameter4', [0 1], 0.001
};

% set initial parameter values in struct objFctParams 
objFctParams.parameter1 =  27.5871;
objFctParams.parameter2 = 88.5443;
objFctParams.parameter3 = 0.7639;
objFctParams.parameter4 = 0.1;

% set single additional function parameter
objFctSettings = 100;

% get default DE parameters
DEParams = getdefaultparams;

% set number of population members (often 10*D is suggested) 
DEParams.NP = 100;

% do not use slave process here
DEParams.feedSlaveProc = 0;

% set times
DEParams.maxiter  = 20;
DEParams.maxtime  = 30; % in seconds
DEParams.maxclock = [];

% set display options
DEParams.infoIterations = 1;
DEParams.infoPeriod     = 10; % in seconds

% do not send E-mails
emailParams = [];

% set random state in order to always use the same population members here
setrandomseed(1);

% % --- Setting current directory ---
% cd(curr_dir)

% start differential evolution
[bestmem, bestval, bestFctParams, nrOfIterations, resultFileName] = differentialevolution(...
	DEParams, paramDefCell, objFctHandle, objFctSettings, objFctParams, emailParams, optimInfo); %#ok

disp(' ');
disp('Best parameter set returned by function differentialevolution:');
disp(bestFctParams);

% continue optimization by loading result file
if DEParams.saveHistory
  
  disp(' ');
  disp(textwrap2(sprintf(...
    'Now continuing optimization by loading result file %s.', resultFileName)));
  disp(' ');
  
  DEParams.maxiter = 100;
  DEParams.maxtime = 60; % in seconds

  [bestmem, bestval, bestFctParams] = differentialevolution(...
    DEParams, paramDefCell, objFctHandle, objFctSettings, objFctParams, emailParams, optimInfo, ...
    resultFileName); %#ok
  
  disp(' ');
  disp('Best parameter set returned by function differentialevolution:');
  disp(bestFctParams);
end

