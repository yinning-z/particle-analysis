function xy_grad                                ... gradations and summary

%
%   Parameters
%

% specimen IDs
nStrt = 1;
nEnd = nStrt;

%%% Sieve series [mm]
%       #3½ #5  #10 #18 #35  #60  #120   #230                               % QF experimental/simplified sieve series
%        |   |   |   |   |    |     |      |
svs = [ 5.6  4   2   1  0.5  0.25  0.125  0.063 0.04  0.02  0.01  0.005  0.0025  0.00125  0.0006  0.0003 ]; % sieve size [mm]

gradSumName = '[grad]'; % summary folder (out)

% _________________________________________________________________________

%
%   Start
%

% number of samples and sieves
nSmps = nEnd-nStrt+1;
nSvs = size(svs,2);

% area retained per sieve, and table with area gradation
aRet = zeros(nSvs+1,nSmps);
grads = zeros(nSvs+1,nSmps);

% summary of shape parameters. 1000 > max(nAgs) in specimens
AIs = NaN(1000,nSmps+1);
FIs = NaN(1000,nSmps+1);
AIs(:,1) = (1:1000)';
FIs(:,1) = (1:1000)';
shape = 0;

% create directory for plots
if ~exist(gradSumName,'dir'), mkdir(gradSumName); end

% specimen-wise
for n = nStrt:nEnd
    
    dirName = num2str(n);
    
    %%% Read and sort table (aggregates)
    
    % read table 
    tableAg = load([dirName '/AG_table.txt'],'-ascii');
    nAgs = size(tableAg,1);
    
    % shape parameters summary
    if size(tableAg,2) == 10
        % first specimen controls
        if n == nStrt, shape = 1; end
        
        % shape parameters have been calculated
        AIs(1:nAgs,n+1) = tableAg(:,9); % AIs
        FIs(1:nAgs,n+1) = tableAg(:,10); % FIs
    end
    
    % sort table by charLen (sieve size)
    [~,index] = sort(tableAg(:,6),'descend');  % Col. 6, charLen
    tableAg = tableAg(index,:);
    %tableAg(:,1)=(1:nAgs)';
    
    %%% Calculate area retained per sieve
    
    % current sieve
    sv = 1;
        
    % polygon-wise
    for i = 1:nAgs
        
        % area and charLen of the polygon
        area = tableAg(i,3);
        charLen = tableAg(i,6);
        
        % find retaining sieve from series (not from table)
        added = 0;
        while ~added
            
            if sv > nSvs || charLen >= svs(sv)
                aRet(sv,n) = aRet(sv,n) + area;
                added = 1;
            else
                sv = sv+1;
            end
            
        end
        
    end
    
    %%% Calculate gradation (percentage passing)
    
    % retained area normalized as percentage of total area
    grads(:,n) = aRet(:,n)/sum(tableAg(:,3))*100;
    
    % sieve-wise, calculate percentage passing
    grads(1,n) = 100-grads(1,n);
    for i = 2:nSvs+1
        grads(i,n) = grads(i-1,n)-grads(i,n);
    end
    
    % display summary
    agArea = sum(tableAg(:,3)); % [mm²]
    disp(['   Specimen ' num2str(n) '. nAgs = ' num2str(nAgs) ', agArea = ' num2str(agArea) ' mm².' ]);

end

%%% Summarize nVrts and nAgs per sieve

% summary of nVrts
nVrtsMu = zeros(nSvs+1,nSmps);
nVrtsSigma = zeros(nSvs+1,nSmps);

% summary of nAgs per sieve
nAgsSvs = zeros(nSvs+1,nSmps);

% specimen-wise
for n = nStrt:nEnd
    
    dirName = num2str(n);
    
    %%% Read and sort table (aggregates)
    
    % read table 
    tableAg = load([dirName '/AG_table.txt'],'-ascii');
    nAgs = size(tableAg,1);
    
    % sort table by charLen (sieve size)
    [~,index] = sort(tableAg(:,6),'descend');  % Col. 6, charLen
    tableAg = tableAg(index,:);
    
    %%% Find sieve changes and update summaries
    
    sv = 1; % current sieve
    agSv = 1; % larger aggregate in current sieve
    
    % polygon-wise
    for i = 1:nAgs
        
        % area and charLen of the polygon
        charLen = tableAg(i,6);
        
        % find smallest retaining sieve from series (not from table)
        ret = 0;
        while ~ret && sv < nSvs
            
            % find change in sieve
            if charLen < svs(sv) && charLen > svs(sv+1)
                % within sieve
                ret = 1;
            else
                % next sieve
                sv = sv+1;
            end
        end
        if i == 1
            svAnt = sv; % previous sieve
        end
        
        if sv > nSvs % PAN
            nVrtsMu(sv,n) = mean(tableAg(i:end,2));
            nVrtsSigma(sv,n) = std(tableAg(i:end,2));
            
            nAgsSvs(sv,n) = nAgs-i+1;
            
            break;
        end
        
        nAgsSvs(sv+1,n) = nAgsSvs(sv+1,n)+1;
        
        % do something
        if svAnt < sv
            nVrtsMu(sv,n) = mean(tableAg(agSv:i-1,2));
            nVrtsSigma(sv,n) = std(tableAg(agSv:i-1,2));
            agSv = i;
            svAnt = sv;
        end
        
        if i == nAgs % last ag
            nVrtsMu(sv+1,n) = mean(tableAg(agSv:i,2));
            nVrtsSigma(sv+1,n) = std(tableAg(agSv:i,2));
        end
        
    end
    
end


%%% Save summary files

if shape
    % save shape parameters summary
    save([gradSumName '/AIs.txt'], 'AIs', '-ascii');
    save([gradSumName '/FIs.txt'], 'FIs', '-ascii');
end

% save gradation summary
aRet = [([svs 0])' aRet];
grads = [([svs 0])' grads];
save([gradSumName '/aRet.txt'], 'aRet', '-ascii');
save([gradSumName '/grads.txt'], 'grads', '-ascii');

% save nVrts and nAgs summary
nVrtsMu = [([svs 0])' nVrtsMu];
nVrtsSigma = [([svs 0])' nVrtsSigma];
save([gradSumName '/nVrtsMu.txt'], 'nVrtsMu', '-ascii');
save([gradSumName '/nVrtsSigma.txt'], 'nVrtsSigma', '-ascii');

nAgsSvs = [([svs 0])' nAgsSvs];
save([gradSumName '/nAgsSvs.txt'], 'nAgsSvs', '-ascii');

end