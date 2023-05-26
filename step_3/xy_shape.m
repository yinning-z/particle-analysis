function xy_shape                                             ... ag AI/FI

%
%   Parameters
%

% specimen IDs
nStrt = 1;
nEnd = nStrt;

% fixed directional angle
dThetaDeg = 4; % [°]                                                        % REF 1

% angularity equation settings
setsAI = 2;
switch setsAI
    
    case 1
        step = 1; % angularity summation step
        % equation: normalizing, units [radians]
        
    case 2 % used by Al-Rousan et al.
        step = 3; % angularity summation step                               % REF 1
        % equation: not normalizing, units [degrees]
end

% progress display info
pcInc = 0.25; % 25-percent increments

% internal plot + save (time-consuming)
pltInt = 0; % [0-1] plot
pausePlt = 0; % [0-1] pause after each polygon
plotDirName = '[plot]'; % name of plot directory

% rotate aggregate to horizontal position before calculating indexes
horizontalize = 1;

% relax criteria to allow calculation for more shapes
relax = 1;

%__________________________________________________________________________

% 
%   Start
% 

% fixed directional angles
dTheta = dThetaDeg*pi/180; % [rad]
nLines = 2*pi/dTheta; % no repetition at 0°. 90 lines for dtheta = 4°

% specimen-wise
for n = nStrt:nEnd
    
    disp(' '); disp([' §   Specimen ' num2str(n)]); disp(' ');

%   _______________________________________________________________________  ...  1  ...
%   Read files and build cells
    disp('     Reading files...');
    
    dirName = num2str(n);
    
    % read AG table and polygon coordinate file
    tableAg = load([dirName '/AG_table.txt'], '-ascii');
    xysAg = load([num2str(n) '/AG_xy.txt'], '-ascii');
    nAgs = max(xysAg(:,1));
    
    % cell containing polygon coordinates
    ags = cell(nAgs,1); 
    
    lnA = 1; % start line, current polygon
    
    % polygon-wise, read coordinates
    for i = 1:nAgs
        nVrt = tableAg(i,2);
        lnB = lnA+nVrt; % end line, current polygon
        
        % read polygon and add to the cell
        xy = xysAg(lnA:lnB, 2:3); % closed
        ags{i} = xy; % closed
        
        lnA = lnB+2; % start line, next polygon
    end
    
    % read centres
    AG_xyc = load([dirName '/AG_xyc.txt'], '-ascii');
    xcs = AG_xyc(:,1);
    ycs = AG_xyc(:,2);
    
    disp(['     Building tables and calculating AI/FI ' ...
          '(' num2str(nAgs) ' aggregates)...']); disp(' ');
    
    % create directory for plots if needed
    if pltInt && ~exist(plotDirName,'dir'), mkdir(plotDirName); end
    
    % polygon-wise
    for i = 1:nAgs
        
%       ___________________________________________________________________  ...  2  ...
%       Recall coordinates and rotate to horizontal position
        
        % recall aggregate coordinates, 'rotated'
        xy = ags{i}; % closed
        x = xy(:,1);
        y = xy(:,2);
        
        % center at zero
        x = x-xcs(i); % centroid
        y = y-ycs(i); % centroid
       %x = x-mean(x(1:end-1));
       %y = y-mean(y(1:end-1));
        
        switch pltInt; case 1 % plot 'rotated' polygon                      %
            figure(1); clf;                                                 %
            plot(x,y,'k:');                                                 %
            xlabel('x [mm]'); ylabel('y [mm]');                             %
            title(['Polygon ' num2str(i) '. '                          ...  %
                   '\Delta\theta = ' num2str(dThetaDeg) '°, '          ...  %
                   'nLines = ' num2str(nLines)]);                           %
            axis equal; hold on;                                            %
        end                                                                 %
        
        switch horizontalize; case 1
            
            % rotate polygon to horizontal position
            theta = tableAg(i,8);
            R = [cos(theta),sin(theta);-sin(theta),cos(theta)]; % clockwise, MG specimens
            %R = [cos(theta),-sin(theta);sin(theta),cos(theta)]; % counter-clockwise
            xyRot = R*[x';y'];
            x = xyRot(1,:)';
            y = xyRot(2,:)';
            
            switch pltInt; case 1 % plot horizontal polygon                 %
                plot(x,y,'b');                                              %
            end                                                             %
            
        end
        
%       ___________________________________________________________________  ...  3  ...
%       Create table with intersect coordinates, radius and angle

        % table of directional angle info (x,y,r,theta)
        xyrt = zeros(nLines+1,4);
        
        % define length of intersecting/rotating line
        rMax = tableAg(i,4);
        canon = 1.5*rMax; % rMax
        
        %%% Fill in table
        
        theta = 0;
        ok = 1;
        
        % angle-wise
        for j = 1:nLines
            
            % test that centroid is inside polygon
            if j == 1
                in = inpolygon(xcs(i),ycs(i), ags{i}(:,1),ags{i}(:,2));
                if ~in
                    ok = 0;
                    switch pltInt, case 1                                   %
                        plot(xcs(i),ycs(i),'kx')                            %
                    end
                    break
                end
            end
            
            % line at angle theta
            xLine = [0, canon*cos(theta)];
            yLine = [0, canon*sin(theta)];
            
            % line intersection with polygon
            [xInt,yInt] = polyxpoly(xLine,yLine,x,y);
            
            switch pltInt; case 1 % plot and number lines                   %
                plot(xLine,yLine,':','color',[.7 .7 .7]);                   %
                text(xLine(2),yLine(2),num2str(j),                     ...  %
                     'FontSize',8, 'HorizontalAlignment','center');         %
            end                                                             %
            
            % check for validity
            if isscalar(xInt)
                
                % unique intersection [ok]
                r = sqrt(xInt^2+yInt^2);
                xyrt(j,:) = [xInt, yInt, r, theta];
                
                switch pltInt; case 1 % plot intersection                   %
                    plot(xInt,yInt,'b.');                                   %
                end
                
            elseif relax
                
                % several intersections [relaxed, ok]
                rkMax = zeros(1,length(xInt));
                for k = 1:length(xInt)
                    rkMax(k) = sqrt(xInt(k)^2+yInt(k)^2);                    
                end
                [r,ind] = max(rkMax);
                
                xyrt(j,:) = [xInt(ind), yInt(ind), r, theta];
                
                switch pltInt; case 1 % plot intersection                   %
                    plot(xInt,yInt,'rx');                                   %
                    plot(xInt(ind),yInt(ind),'bo');
                end
                
            else
                
                % multiple intersections exist [~ok]

                % Possible reasons include:
                % . polygon coordinates may be unevenly distributed on edges (coordinate centre != centre of mass)
                % . polygon may be non-convex, or polygon centre may lay outside polygon

                ok = 0;

                switch pltInt; case 1 % plot intersection (failed)          %
                    plot(xInt,yInt,'rx');                                   %
                end                                                         %

                break
            
            end
            
            % next angle
            theta = theta+dTheta;

        end
        
        % last line with theta equal to one complete rotation
        xyrt(nLines+1,1:3) = xyrt(1,1:3);
        xyrt(nLines+1,4) = 2*pi;
        
%       ___________________________________________________________________  ...  4  ...
%       Calculate shape parameters, update table and save
        
        if ok
            
            %%% Calculate gradient values for AI
            
            % vector of gradient values
            grads = zeros(nLines,1);
            
            % line-wise
            for j = 1:nLines
                
                %%% Calculate inverse tangent and adjust gradient
                
                x1 = xyrt(j,1);
                x2 = xyrt(j+1,1);
                y1 = xyrt(j,2);
                y2 = xyrt(j+1,2);
                val = atan((y2-y1)/(x2-x1));

                % discriminate by direction to calculate appropriate gradient
                if y2 >= y1
                    if x2 >= x1
                        % segment pointing towards 1st quadrant
                        towards = 1;
                        grad = val;
                    else % x2 < x1
                        % segment pointing towards 2d quadrant
                        towards = 2;
                        grad = pi+val; % val < 0
                    end
                else % y2 < y1
                    if x2 < x1
                        % segment pointing towards 3rd quadrant
                        towards = 3;
                        grad = pi+val;
                    else % x2 >= x1
                        % segment pointing towards 4rth quadrant
                        towards = 4;
                        grad = 2*pi+val; % val<0
                    end
                end

                % then discriminate by quadrant to adjust full rotations
                theta = xyrt(j,4);
                if theta > pi && theta <= 3*pi/2
                    % 3rd quadrant
                    if towards == 1
                        grad = grad+2*pi;
                    end
                elseif theta > 3*pi/2
                    % 4th quadrant
                    if towards == 1 || towards == 2
                        grad = grad+2*pi;
                    end
                end

                % simplified calculation                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if x1 <= x2
                    % segment pointing towards 1st or 4th quadrant
                    if theta <= pi/2 % 1st quadrant
                        gradX = val;
                    else % theta > pi/2. 2nd, 3rd or 4th quadrant
                        gradX = 2*pi+val;
                    end
                else % x1 > x2. segment pointing towards 2nd or 3rd quadrant
                    gradX = pi+val;
                    if theta > 3*pi/2 % 4th quadrant
                        gradX = gradX+2*pi;
                    end
                end
                if grad-gradX ~= 0, warning('difference'); end
                %disp([i, j, theta*180/pi, val*180/pi, towards, grad*180/pi, gradX*180/pi]);

                grads(j,1) = grad;
            end

            %%% Calculate angularity index, AI                              % REF 1, EQ 3
            
            angIndx = 0;
            for j = 1:nLines-step
                angIndx = angIndx + abs(grads(j)-grads(j+step));
            end
            
            switch setsAI
                case 1
                    angIndx = angIndx/(nLines/step-1); % normalizing, units [radians]
                case 2
                    angIndx = angIndx*180/pi; % not normalizing, units [degrees]
            end
            
            tableAg(i,9) = angIndx; % AI
            
            %%% Calculate form index, FI                                    % REF 1, EQ 2
            
            formIndx = 0;
            for j = 1:nLines % all lines included. no line repetition at 0°
                rt = xyrt(j,3);
                rtdt = xyrt(j+1,3);
                formIndx = formIndx+abs(rtdt-rt)/rt;
            end
            tableAg(i,10) = formIndx; % FI
        
        else % ~ok
            tableAg(i,9) = NaN; % AI
            tableAg(i,10) = NaN; % FI
            
            disp(['     * Specimen ' num2str(n) ', '                   ...
                  'AI/FI not calculated for polygon ' num2str(i) '.']);
        end
        
%       ___________________________________________________________________  ...  5  ...
%       Display progress/update plot and save, pause if requested
        
        switch pltInt; case 1
            
            % display AI/FI, save                                           %
            if ok                                                           %
                text(0,0, ['AI = ' num2str(angIndx)],                  ...  %
                     'HorizontalAlignment','center', 'VerticalAlignment','bottom');
                text(0,0, ['FI = ' num2str(formIndx)],                 ...  %
                     'HorizontalAlignment','center', 'VerticalAlignment','top');
                saveas(gcf,[plotDirName '/' dirName '-' num2str(i) ', '...  %
                       'AI = ' num2str(angIndx) ', '                   ...  %
                       'FI = ' num2str(formIndx)  '.png']);                 %
            else % ~ok                                                      %
                saveas(gcf,[plotDirName '/[NaN] ' dirName '-' num2str(i) '.png']);
            end                                                             %
            
            % pause if requested
            switch pausePlt
                case 0
                    drawnow limitrate
                case 1
                    if n == 1 && i == 1
                        disp('     Press any key to continue to the next polygon... ');
                        disp('     To stop, click on the Command Window and press Ctrl+C ');
                    end
                    if ~(n == nEnd && i == nAgs)
                        pause;
                    end
            end
            
        otherwise % ~pltInt
            
            % display progress
            if mod(i,round(pcInc*nAgs)) == 0 || i == nAgs
                disp([ num2str(round(100*i/nAgs)) ' %' ]);
            end
            
        end
        
    end
    
    % save updated table
    save([dirName '/AG_table.txt'], 'tableAg', '-ascii');
    
    disp(' ');

end

end

%{
    __________________________________________________________________
    REF

 1. "New Methodology for Shape Classification of Aggregates"
     Al-Rousan T, Masad E, Myers L, Speigelman C [2005].

  * dThetadeg = 4°, from Section "Shape analysis methods", Form.
  * AI (gradient method), modified from Equation 3.
  * FI, Equation 2.

%}