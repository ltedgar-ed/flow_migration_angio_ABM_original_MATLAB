% Plot the vessel network along with some pressure and flow values (panel 1)
% Plot the cell polarity vectors (panel 2)
% subfunction of ABM_EC_polarize_migrate_flow_H.m
% 08/20/2018
% --------------------
function plot_network_Ub1(segments, D, P, Q, tau, seg_cells)

% Close plot if not the first time plotting
if (ishandle(1))
    close(1)
end

%% Pressure, Flow, and Diameter Plot
figure (1), subplot(1,3,1),
set(gcf, 'Position', [60, 290, 1800, 500])

hold on
axis([-50 150 -30 120])
axis square

% Split segments1 into three vessels
segments1 = segments(:,1:2);

vessel1 = segments1(1:6,:);
vessel2 = segments1(7:17,:);
vessel3 = segments1(18:23,:);

% Plot the first vessel
for seg = 1:5
    if (D(seg) ~= 0)
        plot([vessel1(seg,1) vessel1(seg+1,1)], [vessel1(seg,2) vessel1(seg+1,2)], 'r-','MarkerSize',14,'LineWidth',D(seg)*1e6/2)
    end
end

% Plot the second vessel
for seg = 6:15
    if (D(seg) ~= 0)
        plot([vessel2(seg-5,1) vessel2(seg-5+1,1)], [vessel2(seg-5,2) vessel2(seg-5+1,2)], 'r-','MarkerSize',14,'LineWidth',D(seg)*1e6/2)
    end
end

% Plot the third vessel
for seg = 16:20
    if (D(seg) ~= 0)
        plot([vessel3(seg-15,1) vessel3(seg-15+1,1)], [vessel3(seg-15,2) vessel3(seg-15+1,2)], 'r-','MarkerSize',14,'LineWidth',D(seg)*1e6/2)
    end
end

% Split segments2 into three vessels
segments2 = segments(:,3:4);

vessel4 = segments2(1:6,:);
vessel5 = segments2(7:17,:);
vessel6 = segments2(18:23,:);

% Plot the fourth vessel
for seg = 1:5
    if (D(seg+20) ~= 0)
        plot([vessel4(seg,1) vessel4(seg+1,1)], [vessel4(seg,2) vessel4(seg+1,2)], 'r-','MarkerSize',14,'LineWidth',D(seg+20)*1e6/2)
    end
end

% Plot the fifth vessel
for seg = 6:15
    if (D(seg+20) ~= 0)
        plot([vessel5(seg-5,1) vessel5(seg-5+1,1)], [vessel5(seg-5,2) vessel5(seg-5+1,2)], 'r-','MarkerSize',14,'LineWidth',D(seg+20)*1e6/2)
    end
end

% Plot the sixth vessel
for seg = 16:20
    if (D(seg+20) ~= 0)
        plot([vessel6(seg-15,1) vessel6(seg-15+1,1)], [vessel6(seg-15,2) vessel6(seg-15+1,2)], 'r-','MarkerSize',14,'LineWidth',D(seg+20)*1e6/2)
    end
end

% Plot the network backbone
plot(segments1(:,1), segments1(:,2), 'k:','LineWidth',0.1)
plot(segments2(:,1), segments2(:,2), 'k:','LineWidth',0.1)

% Display pressure along the network (in black)
text(-8,-5,['Pin = ' num2str(P(1)/98)],'Color','k','FontSize',10)
text(100-8,-5,['Pout = ' num2str(P(21)/98)],'Color','k','FontSize',10)

text(-25,50,[num2str(P(6)/98)],'Color','k','FontSize',10)
text(47,60,[num2str(P(11)/98)],'Color','k','FontSize',10)
text(105,50,[num2str(P(16)/98)],'Color','k','FontSize',10)

text(-25,100,[num2str(P(26)/98)],'Color','k','FontSize',10)
text(47,110,[num2str(P(31)/98)],'Color','k','FontSize',10)
text(105,100,[num2str(P(36)/98)],'Color','k','FontSize',10)

% Display flow along the network (in blue)
text(-30,25,[num2str(Q(3)*3.6e12)],'Color','b','FontSize',9)
text(-30,75,[num2str(Q(23)*3.6e12)],'Color','b','FontSize',9)
text(105,25,[num2str(Q(18)*3.6e12)],'Color','b','FontSize',9)
text(105,75,[num2str(Q(38)*3.6e12)],'Color','b','FontSize',9)

text(15,55,[num2str(Q(8)*3.6e12)],'Color','b','FontSize',9)
text(65,55,[num2str(Q(13)*3.6e12)],'Color','b','FontSize',9)

text(15,105,[num2str(Q(28)*3.6e12)],'Color','b','FontSize',9)
text(65,105,[num2str(Q(33)*3.6e12)],'Color','b','FontSize',9)

% shear stress
text(5,25,[num2str(tau(3))],'Color','r','FontSize',9)
text(5,75,[num2str(tau(23))],'Color','r','FontSize',9)
text(70,25,[num2str(tau(18))],'Color','r','FontSize',9)
text(70,75,[num2str(tau(38))],'Color','r','FontSize',9)

text(15,45,[num2str(tau(8))],'Color','r','FontSize',9)
text(65,45,[num2str(tau(13))],'Color','r','FontSize',9)

text(15,95,[num2str(tau(28))],'Color','r','FontSize',9)
text(65,95,[num2str(tau(33))],'Color','r','FontSize',9)

text(-45,-25,['Pressure, P (cmH2O)'],'Color','k','FontSize',10)
text(45,-25,['Flow, Q (\muL/hr)'],'Color','b','FontSize',10)
text(110,-25,['WSS (Pa)'],'Color','r','FontSize',10)

grid on
title(' Pressure, Flow, Diameter of Network ')
xlabel(' X (\mum) ')
ylabel(' Y (\mum) ')
hold off

%% Cell Polarity Plot
% Plot the second panel
figure (1), subplot(1,3,2),
hold on
axis([-50 150 -30 120])
axis square

cell_space = 2;                                                             % Space between cell nuclei
S = 12;                                                                     % Size of cell nuclei

L = 4;                                                                      % Length of polarity vector
W = 0.75;                                                                   % Width of polarity vector
H = 8;                                                                      % Size of Golgi apparatus 

% Plot the first vessel
for seg = 1:5
    if (seg_cells{seg,1} ~= 0)
        x_span = (seg_cells{seg,1} - 1)*cell_space;
        mid_x = (vessel1(seg,1) + vessel1(seg+1,1))/2;
        mid_y = (vessel1(seg,2) + vessel1(seg+1,2))/2;
        
        % If only one cell in the segment...
        if (seg_cells{seg,1} == 1)
            plot(mid_x, mid_y, 'r.','MarkerSize',S)                         % Plot cell nuclei
            
            u = seg_cells{seg,2}(1,1);
            v = seg_cells{seg,2}(2,1); 
            plot([(mid_x - x_span/2 - (L/2)*u) (mid_x + (L/2)*u)], [mid_y - (L/2)*v mid_y + (L/2)*v],'b','LineWidth',W)     % Plot the polarity vector
            plot((mid_x + (L/2)*u), mid_y + (L/2)*v, 'b.','MarkerSize',H)   % Plot the Golgi apparatus
        else              
            for cell = 1:seg_cells{seg,1}
                plot(mid_x + (cell-1)*cell_space - x_span/2, mid_y, 'r.','MarkerSize',S)                         % Plot cell nuclei
                
                u = seg_cells{seg,2}(1,cell);
                v = seg_cells{seg,2}(2,cell);               
                plot([(mid_x + (cell-1)*cell_space - x_span/2 - (L/2)*u) (mid_x + (cell-1)*cell_space - x_span/2 + (L/2)*u)], [mid_y - (L/2)*v mid_y + (L/2)*v],'b','LineWidth',W)     % Plot the polarity vector 
                plot((mid_x + (cell-1)*cell_space - x_span/2 + (L/2)*u), mid_y + (L/2)*v, 'b.','MarkerSize',H)   % Plot the Golgi apparatus
            end
        end
    end
end

% Plot the second vessel
for seg = 6:15
    if (seg_cells{seg,1} ~= 0)
        y_span = (seg_cells{seg,1} - 1)*cell_space;
        mid_x = (vessel2(seg-5,1) + vessel2(seg-5+1,1))/2;
        mid_y = (vessel2(seg-5,2) + vessel2(seg-5+1,2))/2;
        
        % If only one cell in the segment...
        if (seg_cells{seg,1} == 1)
            plot(mid_x, mid_y, 'r.','MarkerSize',S)                         % Plot cell nuclei
            
            u = seg_cells{seg,2}(1,1);
            v = seg_cells{seg,2}(2,1); 
            plot([mid_x - (L/2)*u mid_x + (L/2)*u], [(mid_y - y_span/2 - (L/2)*v) (mid_y + (L/2)*v)] ,'b','LineWidth',W)     % Plot the polarity vector 
            plot(mid_x + (L/2)*u, (mid_y + (L/2)*v),  'b.','MarkerSize',H)   % Plot the Golgi apparatus
        else              
            for cell = 1:seg_cells{seg,1}
                plot(mid_x, mid_y + (cell-1)*cell_space - y_span/2, 'r.','MarkerSize',S)                         % Plot cell nuclei
                
                u = seg_cells{seg,2}(1,cell);
                v = seg_cells{seg,2}(2,cell);               
                plot([mid_x - (L/2)*u mid_x + (L/2)*u], [(mid_y + (cell-1)*cell_space - y_span/2 - (L/2)*v) (mid_y + (cell-1)*cell_space - y_span/2 + (L/2)*v)],'b','LineWidth',W)     % Plot the polarity vector 
                plot(mid_x + (L/2)*u, (mid_y + (cell-1)*cell_space - y_span/2 + (L/2)*v), 'b.','MarkerSize',H)   % Plot the Golgi apparatus
            end
        end
    end
end

% Plot the third vessel
for seg = 16:20
    if (seg_cells{seg,1} ~= 0)
        x_span = (seg_cells{seg,1} - 1)*cell_space;
        mid_x = (vessel3(seg-15,1) + vessel3(seg-15+1,1))/2;
        mid_y = (vessel3(seg-15,2) + vessel3(seg-15+1,2))/2;
        
        % If only one cell in the segment...
        if (seg_cells{seg,1} == 1)
            plot(mid_x, mid_y, 'r.','MarkerSize',S)                         % Plot cell nuclei
            
            u = seg_cells{seg,2}(1,1);
            v = seg_cells{seg,2}(2,1); 
            plot([(mid_x - x_span/2 - (L/2)*u) (mid_x + (L/2)*u)], [mid_y - (L/2)*v mid_y + (L/2)*v],'b','LineWidth',W)     % Plot the polarity vector 
            plot((mid_x + (L/2)*u), mid_y + (L/2)*v, 'b.','MarkerSize',H)   % Plot the Golgi apparatus
        else              
            for cell = 1:seg_cells{seg,1}
                plot(mid_x + (cell-1)*cell_space - x_span/2, mid_y, 'r.','MarkerSize',S)                         % Plot cell nuclei
                
                u = seg_cells{seg,2}(1,cell);
                v = seg_cells{seg,2}(2,cell);               
                plot([(mid_x + (cell-1)*cell_space - x_span/2 - (L/2)*u) (mid_x + (cell-1)*cell_space - x_span/2 + (L/2)*u)], [mid_y - (L/2)*v mid_y + (L/2)*v],'b','LineWidth',W)     % Plot the polarity vector 
                plot((mid_x + (cell-1)*cell_space - x_span/2 + (L/2)*u), mid_y + (L/2)*v, 'b.','MarkerSize',H)   % Plot the Golgi apparatus
            end
        end
    end
end

% Plot the fourth vessel
for seg = 21:25
    if (seg_cells{seg,1} ~= 0)
        x_span = (seg_cells{seg,1} - 1)*cell_space;
        mid_x = (vessel4(seg-20,1) + vessel4(seg-20+1,1))/2;
        mid_y = (vessel4(seg-20,2) + vessel4(seg-20+1,2))/2;
        
        % If only one cell in the segment...
        if (seg_cells{seg,1} == 1)
            plot(mid_x, mid_y, 'r.','MarkerSize',S)                         % Plot cell nuclei
            
            u = seg_cells{seg,2}(1,1);
            v = seg_cells{seg,2}(2,1); 
            plot([(mid_x - x_span/2 - (L/2)*u) (mid_x + (L/2)*u)], [mid_y - (L/2)*v mid_y + (L/2)*v],'b','LineWidth',W)     % Plot the polarity vector
            plot((mid_x + (L/2)*u), mid_y + (L/2)*v, 'b.','MarkerSize',H)   % Plot the Golgi apparatus
        else              
            for cell = 1:seg_cells{seg,1}
                plot(mid_x + (cell-1)*cell_space - x_span/2, mid_y, 'r.','MarkerSize',S)                         % Plot cell nuclei
                
                u = seg_cells{seg,2}(1,cell);
                v = seg_cells{seg,2}(2,cell);               
                plot([(mid_x + (cell-1)*cell_space - x_span/2 - (L/2)*u) (mid_x + (cell-1)*cell_space - x_span/2 + (L/2)*u)], [mid_y - (L/2)*v mid_y + (L/2)*v],'b','LineWidth',W)     % Plot the polarity vector 
                plot((mid_x + (cell-1)*cell_space - x_span/2 + (L/2)*u), mid_y + (L/2)*v, 'b.','MarkerSize',H)   % Plot the Golgi apparatus
            end
        end
    end
end

% Plot the fifth vessel
for seg = 26:35
    if (seg_cells{seg,1} ~= 0)
        y_span = (seg_cells{seg,1} - 1)*cell_space;
        mid_x = (vessel5(seg-25,1) + vessel5(seg-25+1,1))/2;
        mid_y = (vessel5(seg-25,2) + vessel5(seg-25+1,2))/2;
        
        % If only one cell in the segment...
        if (seg_cells{seg,1} == 1)
            plot(mid_x, mid_y, 'r.','MarkerSize',S)                         % Plot cell nuclei
            
            u = seg_cells{seg,2}(1,1);
            v = seg_cells{seg,2}(2,1); 
            plot([mid_x - (L/2)*u mid_x + (L/2)*u], [(mid_y - y_span/2 - (L/2)*v) (mid_y + (L/2)*v)] ,'b','LineWidth',W)     % Plot the polarity vector 
            plot(mid_x + (L/2)*u, (mid_y + (L/2)*v),  'b.','MarkerSize',H)   % Plot the Golgi apparatus
        else              
            for cell = 1:seg_cells{seg,1}
                plot(mid_x, mid_y + (cell-1)*cell_space - y_span/2, 'r.','MarkerSize',S)                         % Plot cell nuclei
                
                u = seg_cells{seg,2}(1,cell);
                v = seg_cells{seg,2}(2,cell);               
                plot([mid_x - (L/2)*u mid_x + (L/2)*u], [(mid_y + (cell-1)*cell_space - y_span/2 - (L/2)*v) (mid_y + (cell-1)*cell_space - y_span/2 + (L/2)*v)],'b','LineWidth',W)     % Plot the polarity vector 
                plot(mid_x + (L/2)*u, (mid_y + (cell-1)*cell_space - y_span/2 + (L/2)*v), 'b.','MarkerSize',H)   % Plot the Golgi apparatus
            end
        end
    end
end

% Plot the sixth vessel
for seg = 36:40
    if (seg_cells{seg,1} ~= 0)
        x_span = (seg_cells{seg,1} - 1)*cell_space;
        mid_x = (vessel6(seg-35,1) + vessel6(seg-35+1,1))/2;
        mid_y = (vessel6(seg-35,2) + vessel6(seg-35+1,2))/2;
        
        % If only one cell in the segment...
        if (seg_cells{seg,1} == 1)
            plot(mid_x, mid_y, 'r.','MarkerSize',S)                         % Plot cell nuclei
            
            u = seg_cells{seg,2}(1,1);
            v = seg_cells{seg,2}(2,1); 
            plot([(mid_x - x_span/2 - (L/2)*u) (mid_x + (L/2)*u)], [mid_y - (L/2)*v mid_y + (L/2)*v],'b','LineWidth',W)     % Plot the polarity vector 
            plot((mid_x + (L/2)*u), mid_y + (L/2)*v, 'b.','MarkerSize',H)   % Plot the Golgi apparatus
        else              
            for cell = 1:seg_cells{seg,1}
                plot(mid_x + (cell-1)*cell_space - x_span/2, mid_y, 'r.','MarkerSize',S)                         % Plot cell nuclei
                
                u = seg_cells{seg,2}(1,cell);
                v = seg_cells{seg,2}(2,cell);               
                plot([(mid_x + (cell-1)*cell_space - x_span/2 - (L/2)*u) (mid_x + (cell-1)*cell_space - x_span/2 + (L/2)*u)], [mid_y - (L/2)*v mid_y + (L/2)*v],'b','LineWidth',W)     % Plot the polarity vector 
                plot((mid_x + (cell-1)*cell_space - x_span/2 + (L/2)*u), mid_y + (L/2)*v, 'b.','MarkerSize',H)   % Plot the Golgi apparatus
            end
        end
    end
end

% Plot the network backbone
plot(segments1(:,1), segments1(:,2), 'k:','LineWidth',0.1)
plot(segments2(:,1), segments2(:,2), 'k:','LineWidth',0.1)

grid on
title(' Cell Polarity ')
xlabel(' X (\mum) ')
ylabel(' Y (\mum) ')
hold off

%% Cell Polarity Distribution Plot
% Plot the second panel
figure (1), subplot(1,3,3),
hold on
axis([-1 1 -1 1])
axis square 

% Plot the polarity unit vector for each cell (unit circle for a perfectly random distribution
for seg = 1:30
    for cell = 1:seg_cells{seg,1}
        plot([0 seg_cells{seg,2}(1,cell)], [0 seg_cells{seg,2}(2,cell)], 'b-')
    end
end
 
grid on
title(' Distribution of Cell Polarity ')
xlabel(' u ')
ylabel(' v ')
hold off    
