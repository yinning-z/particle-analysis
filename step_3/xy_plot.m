function xy_plot

% 
%  Parameters
% 

% specimen IDs
nStrt = 1;
nEnd = nStrt;

% plot individual aggregates
plotAgs = [];
% plotAgs = [1,2,3,4,5];
% plotAgs = [53,100,132,140,152,200,245,274,281]; % A
% plotAgs = [79,122,206,215,226,231,283,309,330,353,355]; % B
% plotAgs = [8,22,36,72,75,77,81,94,114,122,156,163,175,185,197,206,209,210,233,255,260,268,279,281,286,314,322,345,351,353,370,427,436,465,495,500,513,549,550,586,590,624,674,690,708,709,768,778,804,831,836,848,849,881,932,933,976,1001,1030,1120,1305,1371,1374,1382,1413,1465,1535,1561]; % C
% plotAgs = [80,106,141,182,204,210,225,234,264,266,378,401,411,486,548,551,583,621,689,696,719,773,846,934,987,998,1161,1181,1239,1302,1365,1542,1551,1708]; % D
% plotAgs = [27,86,98,100,139,152,173,174,176,190,195,214,217,228,233,234,241,281,286,338,344,347,359,379,403,421,430,435,546,552,556,579,610,613,654,701,742,796,999,1166]; % E
% plotAgs = [34,43,51,60,106,134,139,166,197,251,262,340,402,444,502,511,520,529,553,583,602,616,703,740,829,866,869,935,948,951,1010,1115,1152,1165,1227,1263,1296,1308,1359,1509,1617,1666]; % F
% plotAgs = [39,44,103,108,123,135,148,190,235,276,306,313,318,321,352,367,406,421,436,460,479,487,534,538,567,635,677,686,699,760,795,808,894,923,970,1000,1028,1119,1151,1199,1895,1982]; % H

% title and axes
tax = 1; % Y/N

% plot directory (output)
dirPlot = '[plot]';

% figure number 
nFigSp = 1; % specimen plot
nFigPols = 2; % polygons plot

% _________________________________________________________________________
% 
%  Start
% 

% create plot directory
if ~exist(dirPlot,'dir'), mkdir(dirPlot); end

% specimen-wise
for n = nStrt:nEnd
    
%   _______________________________________________________________________
%   Read files and build cells

    dirName = num2str(n);
    
    %%% Specimen
    
    % read specimen coordinates
    xySp = load([dirName '/sp_xy.txt'],'-ascii');
    sx = xySp(:,1);
    sy = xySp(:,2);
    
    spArea = polyarea(sx,sy); % specimen area [mm²]
    
    
    %%% Aggregates
    
    % read table
    tableAg = load([dirName '/AG_table.txt'], '-ascii');
    nAgs = size(tableAg,1);
    
    % read polygon coordinate file, update number of polygons
    xysAg = load([num2str(n) '/AG_xy.txt'], '-ascii');
        
    % cell containing polygon coordinates
    ags = cell(nAgs,1);
    
    % polygon-wise
    lnA = 1; % current polygon start line
    for i = 1:nAgs
        nVrt = tableAg(i,2);
        lnB = lnA+nVrt; % current polygon end line
        
        % read polygon and add to the cell
        xy = xysAg(lnA:lnB, 2:3); % closed
        ags{i} = xy; % closed
        
        lnA = lnB+2; % next polygon start line
    end
    
    % read polygon centre coordinates
    cxyAgs = load([dirName '/AG_xyc.txt']);
    cxAg = cxyAgs(:,1);
    cyAg = cxyAgs(:,2);
    
    % display summary
    disp(['   Specimen ' num2str(n) '. nAgs = '  num2str(nAgs) ', agFrac = ' num2str(sum(tableAg(:,3))/spArea*100) '%' ]);
    
%   _______________________________________________________________________
%   Plot

    %%% Plot specimen and aggregates
    
    % figure
    if n == nStrt
        figSp = figure(nFigSp);
    end
    set(0, 'CurrentFigure',figSp)
    clf; hold on
    
    % plot specimen
    fill(sx,sy, [.9 .9 .9], 'EdgeColor',[.75 .75 .75]); % concrete
   %fill(sx,sy, [0.025 0.025 0.025], 'EdgeColor','None'); % HMA

    % polygon-wise, plot aggregates
    for i = 1:nAgs
        
        xy = ags{i}; % closed
        x = xy(:,1);
        y = xy(:,2);
        
        if ismember(i,plotAgs)
            
            %%% Highlight specified aggregate
           %fill(x,y,'r', 'EdgeColor','None');
            fill(x,y,[0.8500 0.3250 0.0980], 'EdgeColor','None');
            text(cxAg(i),cyAg(i), num2str(i), 'FontSize', 12, 'Color','w', 'HorizontalAlignment','left');
            
        else
            
            %%% Flat color
            colorVct = [.7 .7 .7];

            %%% Plot aggregate
            fill(x,y,colorVct, 'EdgeColor',[.5 .5 .5]);
           %fill(x,y,colorVct, 'EdgeColor','None');

        end
        
    end
    
    %%% Finalize/save specimen plot
    
    xlabel('x [mm]'); ylabel('y [mm]')
    axis equal tight
    if ~tax
        axis off
    end
    hold off

    if tax
        title(['Specimen ' num2str(n) '. nAgs = ' num2str(nAgs) ', agFrac = ' num2str(100*sum(tableAg(:,3))/spArea) '%, spArea = ' num2str(spArea) 'mm^2'])
    end
    
    if isempty(plotAgs)
        saveas(figSp,[dirPlot '/' num2str(n) ' - ' num2str(nAgs) ' ags, agFrac = ' num2str(100*sum(tableAg(:,3))/spArea) '%, spArea = ' num2str(spArea) 'mm^2' '.png'])
        % saveas(figSp,[dirPlot '/' num2str(n) ' - ' num2str(nAgs) ' ags, agFrac = ' num2str(100*sum(tableAg(:,3))/spArea) '%, spArea = ' num2str(spArea) 'mm^2' '.emf'])
        % saveas(figSp,[dirPlot '/' num2str(n) ' - ' num2str(nAgs) ' ags, agFrac = ' num2str(100*sum(tableAg(:,3))/spArea) '%, spArea = ' num2str(spArea) 'mm^2' '.emf'], 'meta') % requires figSp.Renderer = 'painters';       
    else
        saveas(figSp,[dirPlot '/' num2str(n) ' - ' num2str(nAgs) ' ags, agFrac = ' num2str(100*sum(tableAg(:,3))/spArea) '%, spArea = ' num2str(spArea) 'mm^2' ' (ags)' '.png'])
    end
    
    
    %%% Plot selected aggregates
    
    if size(plotAgs) ~= 0
        % individual aggregates specified
        
        if n == nStrt
            figPols = figure(nFigPols);
            tiledlayout('flow')
        end
        set(0, 'CurrentFigure',figPols); clf;
        
        % polygon-wise among specified polygons
        nIndPols = length(plotAgs);
        for i = 1:nIndPols
            
            agID = plotAgs(i);
            
            xy = ags{agID}; % closed
            x = xy(:,1) - cxAg(agID); % centered
            y = xy(:,2) - cyAg(agID); % centered
            
            % rotate polygon to horizontal position
            theta = tableAg(agID,8);
            R = [cos(theta),sin(theta);-sin(theta),cos(theta)]; % clockwise, MG specimens
            xyRot = R*[x';y'];
            x = xyRot(1,:)';
            y = xyRot(2,:)';
            
            nexttile
            hold on
           %subplot(1,nIndPols,i); % horizontal
           %subplot(nIndPols,1,i); % vertical
           
            % plot polygon, mark centre
            fill(x,y,colorVct, 'EdgeColor',[.5 .5 .5])
           %plot(x,y,'k');
            plot(0,0, 'kx')
            text(0,0, ['  ' num2str(agID)], 'HorizontalAlignment','left', 'FontSize',8)
            
            % plot retaining sieve to the right
           %svRet = tableAg(agID,7);
           %plot([max(x),max(x)], [-svRet/2,svRet/2], 'r--.')
            
            % finalize plot
            xlabel('x [mm]'); ylabel('y [mm]')
           %title(['Specimen ' num2str(n) '. agID = ' num2str(agID)])
           
            axis equal tight
            
           %xlim([-3 3]) % A, B
            xlim([-1.5 1.5]) % C, D, E, F, H
           %ylim([-1.5 1.5])
           %axis square
            
            hold off
        
        end
        
        % save plot
        lblPltAgs = char(join(string(plotAgs),','));
        
        if length(lblPltAgs)>200
            lblPltAgs = [lblPltAgs(1:200) '...'];
        end
        
        saveas(figPols,[dirPlot '/' num2str(n) ' - ags (' lblPltAgs ').png'])
        
    end

end   

end