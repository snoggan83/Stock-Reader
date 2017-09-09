function sh = leadlag(Close,time,n,m,annualscaling,tag)
% -------------------------------------------------------------------------
%
% Calculation of EMA lead/lag ans Sharpe Ratio
%
% Author: John Nilsson, 2016-06-28
% -------------------------------------------------------------------------

% --- Calculate EMA ---
[lead,lag] = movavg(Close,n,m,'e');         % EMA (lead/lag)

s = zeros(size(Close));
s(lead<lag) = -1;
s(lag<lead) = 1;

r = [0 ; s(1:end-1).*diff(Close)];
w = cumsum(r);                              % Cumulative revenue

sh = annualscaling*sharpe(r,0);             % Annual Sharpe ratio

if strcmp(tag,'plot')
    
    % --- Executing plotting ---
    ax(1) = subplot(2,1,1);
    plot(time,[Close lead lag]) 
    title(strcat('LeadLag EMA Results, Annual Sharpe Ratio =',{' '},num2str(sh)))  
    grid on
    
    legend([{'Close'} strcat('Lead =',{' '},num2str(n)) strcat('Lag =',{' '},num2str(m))]) 
    dynamicDateTicks

    ax(2) = subplot(2,1,2);
    plot(time,[s w])
    title(strcat('Final Return =',{' '},num2str(w(end)),{' '},'(',num2str(round(w(end)/Close(1)*100)),'%)')); 
    grid on
    dynamicDateTicks

    
    
    
    linkaxes(ax,'x');
end