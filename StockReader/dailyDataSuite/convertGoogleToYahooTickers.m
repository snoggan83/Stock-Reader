function yahooTickers = convertGoogleToYahooTickers( googleTickers )
% Converts ticker symbols formatted for Yahoo! Finance to those used by
% Google Finance

% Allow for a single character string as an input.  Set flag so that the
% output is a single char string, too.
if ischar(googleTickers)
    isInputChar = true;
    googleTickers = {googleTickers};
else
    isInputChar = false;
end

yahooTickers = cellfun(@convertGToY, googleTickers, 'UniformOutput', false);

if isInputChar
    yahooTickers = yahooTickers{1};
end

end

function yTicker = convertGToY(gTicker)

sepLoc = strfind(gTicker, ':');

if ~isempty(sepLoc)
    gExchange = gTicker(1:sepLoc-1);
    ticker = gTicker(sepLoc+1:end);
    
    % The full lists are at:
    % http://www.google.com/googlefinance/disclaimer/
    % http://finance.yahoo.com/exchanges
    
    % Format: {Yahoo, Google; ...}
    lookup = {...
        'PA',   'EPA';...
        'MI',   'BIT';...
        'AS',   'AMS';...
        'DE',   'ETR';...
        'BR',   'EBR';...
        'L',    'LON';...
        'none', 'NYSE';...
        'none', 'NASDAQ';...
        'OB',   'OTC';...
        'PK',   'PINK';...
        'TO',   'TSE';...
        'V',    'CVE';...
        'ST',   'STO';...
        'CO',   'CPH';...
        'BO',   'BOM';...
        'NS',   'NSE';...
        'SS',   'SHA';...
        'SZ',   'SHE';...
        'TW',   'TPE';...
        'HK',   'HKG';...
        'SI',   'SGX';...
        'KS',   'KRX';...
        'KQ',   'KOSDAQ';...
        'TA',   'TLV';...
        'AX',   'ASX';...
        'NZ',   'NZE';...
        'none', 'MUTF'};
    
    row = strcmp(gExchange, lookup(:,2));
    if ~any(row)
        % We didn't find a match.  Likely an exchange that only Google
        % supports-- return a blank Yahoo! ticker rather than erroring out.
        yTicker = '';
    else
        yExchange = lookup{row, 1};
                
        if strcmp(yExchange, 'none')
            % Some Google tickers can be written like 'NASDAQ:MSFT' while
            % Yahoo! will only accept 'MSFT' (i.e.: no exchange).
            yTicker = ticker;
        else
            yTicker = [ticker '.' yExchange];
        end
    end
    
else
    % No separator: a US-based stock for which the Yahoo! ID is the same as
    % the Google one (hopefully).
    yTicker = gTicker;
    
end
end