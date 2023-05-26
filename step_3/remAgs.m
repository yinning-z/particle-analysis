function remAgs

% 
%  Parameters
% 

smpName = 'A';
% smpName = 'B';

% smpName = 'C';
% smpName = 'D';
% smpName = 'E';
% smpName = 'F';
% smpName = 'H';

% list of aggregates to be removed
switch smpName
    case {'A'}
        dirName = '1';
        remAgs = [53,100,132,140,152,200,245,274,281];
    case {'B'}
        dirName = '2';
        remAgs = [79,122,206,215,226,231,283,309,330,353,355];
        
    case {'C'}
        dirName = '3';
        remAgs = [8,22,36,72,75,77,81,94,114,122,156,163,175,185,197,206,209,210,233,255,260,268,279,281,286,314,322,345,351,353,370,427,436,465,495,500,513,549,550,586,590,624,674,690,708,709,768,778,804,831,836,848,849,881,932,933,976,1001,1030,1120,1305,1371,1374,1382,1413,1465,1535,1561];
    case {'D'}
        dirName = '4';
        remAgs = [80,106,141,182,204,210,225,234,264,266,378,401,411,486,548,551,583,621,689,696,719,773,846,934,987,998,1161,1181,1239,1302,1365,1542,1551,1708];
    case {'E'}
        dirName = '5';
        remAgs = [27,86,98,100,139,152,173,174,176,190,195,214,217,228,233,234,241,281,286,338,344,347,359,379,403,421,430,435,546,552,556,579,610,613,654,701,742,796,999,1166];
    case {'F'}
        dirName = '6';
        remAgs = [34,43,51,60,106,134,139,166,197,251,262,340,402,444,502,511,520,529,553,583,602,616,703,740,829,866,869,935,948,951,1010,1115,1152,1165,1227,1263,1296,1308,1359,1509,1617,1666];
    case {'H'}
        dirName = '8';
        remAgs = [39,44,103,108,123,135,148,190,235,276,306,313,318,321,352,367,406,421,436,460,479,487,534,538,567,635,677,686,699,760,795,808,894,923,970,1000,1028,1119,1151,1199,1895,1982];
        
    otherwise
        remAgs = [];
end


% 
%  Start
% 

% original aggregates
xysAg = load([dirName '/AG_xy.txt'], '-ascii');
nLines = size(xysAg,1);

% vector of aggregates after removal
xysAg_new = zeros(nLines,3); % trimmed later

% check lines and add those that are not in the 'remove' array
lnNew = 1;
for i = 1:nLines
    nAgi = xysAg(i,1);
    if ismember(nAgi,remAgs)
        % row is not added
    elseif i>1 && ismember(xysAg(i-1,1),remAgs)
        % row is not added
    else
        % row is added
        xysAg_new(lnNew,:) = xysAg(i,:);
        lnNew = lnNew+1;
    end
end

% trim
xysAg_new = xysAg_new(1:lnNew-1,:);

% save new polygon coordinate list file
save(['AG_xy ' smpName ' new.txt'],'xysAg_new','-ascii');

end