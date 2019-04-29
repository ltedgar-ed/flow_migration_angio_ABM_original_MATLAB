% A simple agent-based model of endothelial cells polarizing and migrating 
% within a bifurating vessel in response to blood flow
% U branch 1 network model
% Version 2
% --------------------
% Lowell Taylor Edgar
% Usher Institute
% University of Edinburgh
% 4/10/2019
% --------------------


clc
clear all
close all

rng(123456789);
%rng(987654321);

branch_rule = 6;
branch_alpha = 1.0;

file_tag = ['_2019_04_14_new_branch_rule_', num2str(branch_rule), '_ranseed1_alpha_', num2str(branch_alpha)];

% Input parameters
Nt = 40;                                                                     % Number of time steps
Pin = 4*98;                                                                % Inlet pressure (Pa)
Pout = 1*98;                                                               % Outlet pressure (Pa)

mu = 3.5*1e-3;                                                              % Dynamic viscosity of blood (Pa-s)
    
Nn = 40;                                                                    % Number of nodes
Nseg = 40;                                                                  % Number of segments 

num_cell = 10;                                                               % Initial number of cells per segment

cell_size = 5*1e-6;                                                        % Set the size of each cell (m)

% Polarization re-alignment weights: w1 + w2 + w3 = 1
w2 = 1;                                                                  % Polarization weight - flow component
w3 = 0.00;                                                                  % Polarization weight - re-alignment from neighbors
w4 = 0.00;
w1 = 1 - w2 - w3 - w4;                                                          % Polarization weight - persistence component

% Video file information
file_name = ['ABM_EC_polarize_migrate_flow_Ub1_Ncell', num2str(num_cell), '_Pin_', num2str(Pin/98), '_Pout_', num2str(Pout/98), '_w2_', num2str(w2), '_w3_', num2str(w3),'_w4_', num2str(w4), file_tag];                  % Prepare the video file

vidObj = VideoWriter(file_name,'MPEG-4');
vidObj.FrameRate = 4;
vidObj.Quality = 100;
open(vidObj);

% Initialize arrays
P = zeros(Nn,1);                                                           % Initialize the nodal pressure array (Pa)
Q = zeros(Nseg,1);                                                          % Initialize the segment flow array (m3/s)
G = zeros(Nseg,1);                                                          % Initialize the segment conductance array (m4/Pa-s-m)

H = zeros(Nseg,1);
tau = zeros(Nseg,1);

L = 10*1e-6*ones(Nseg,1);                                                   % Set the segment lengths (m)
Ncell = num_cell*ones(Nseg,1);                                              % Initialize the segment cell number array

segments = make_segments_Ub1(L);                                              % Make the nodal array for plotting

seg_cells = cell(Nseg,3);                                                   % Create the segment cell array
                                                                            %    seg_cells{:,1} - number of cells in the segment
                                                                            %    seg_cells{:,2} - cell polarity vectors
                                                                            %    seg_cells{:,3} - cell migration indicator
% Initialize the segment cell array
for seg = 1:Nseg        
    seg_cells{seg,1} = Ncell(seg);                                          % Number of cells in the segment                                       
    
    for cell = 1:Ncell(seg)
        polar_vect = randn(2,1);                                            % Create a random unit vector for each cell to represent polarity 
        polar_vect = polar_vect/sqrt(polar_vect(1,1)^2 + polar_vect(2,1)^2);
        seg_cells{seg,2} = [seg_cells{seg,2} polar_vect];
    end
    
    seg_cells{seg,3} = zeros(1,Ncell(seg));                                 % Migration indicator, 0 if not migrating
end

% Calculate initial pressure and flow
%D = Ncell*cell_size;                                                        % Initialize the segment diameter array (m)
D = zeros(Nseg,1);

for seg = 1:Nseg
    %if (Ncell(seg) > 1)
    if (Ncell(seg) >= 1)
        D(seg) = Ncell(seg)*cell_size/pi;
    else
        D(seg) = 0;
    end
end

G = (pi*D.^4)./(128*mu*L);                                                  % Calculate segment conductance (m4/Pa-s-m)

for seg = 1:Nseg
    if (D(seg) ~= 0)
        H(seg) = (32*mu)/(pi*D(seg)^3);
    else
        H(seg) = 0;
    end
end

% Solve for flow in the network
[P, Q, tau] = solve_for_flow_Ub1(G, Pin, Pout, H);                     % Solve for nodal pressure (Pa) and segment flow (m3/s)

new_seg_cells = seg_cells;

for seg = 1:Nseg
    [seg_cells, new_seg_cells] = realign_polarity_Ub1(seg, Q, seg_cells, new_seg_cells, w1, w2, w3, w4);
end

% Plot the network
plot_network_Ub1(segments, D, P, Q, tau, seg_cells)                             

% Write to the video file
currFrame = getframe(gcf);
writeVideo(vidObj,currFrame)
  
% Wait for user to click on figure to continue
%waitforbuttonpress                                                          
pause(0.01)

% Begin migration time steps
for t = 1:Nt   
    disp(['Time step ', num2str(t), '/', num2str(Nt)])
    
    migrate = zeros(Nseg,1);                                                % Initialize the segment migratory cell counter
    new_seg_cells = seg_cells;                                              % Create a copy of the segment cell array
    
    % For each segment...
    for seg = 1:Nseg                                                
        
        % Align polarity vectors 
        [seg_cells, new_seg_cells] = realign_polarity_Ub1(seg, Q, seg_cells, new_seg_cells, w1, w2, w3, w4);
        
        [seg_cells, new_seg_cells] = cell_migration_Ub1(seg, seg_cells, new_seg_cells, migrate, Q, tau, branch_rule, branch_alpha);
        
        % Remove old placeholders in polarity vectors array for cells that moved
        new_polar_vects = [];
        [m n] = size(new_seg_cells{seg,2});
        
        for cell = 1:n
            if (norm(new_seg_cells{seg,2}(:,cell)) ~= 0)
                new_polar_vects = [new_polar_vects new_seg_cells{seg,2}(:,cell)];       
            end
        end
        
        new_seg_cells{seg,2} = new_polar_vects;
    
    end
    
    % Update the segment cell array
    seg_cells = new_seg_cells;                                              % Save the new copy of the segment cell array
    
    for seg = 1:Nseg
        Ncell(seg) = seg_cells{seg,1};                                      % Update the segment cell number array
        seg_cells{seg,3} = zeros(1,Ncell(seg));                             % Update the cell migration indicators
    end
    
    % Calculate updated pressure and flow
    %D = Ncell*cell_size;                                                    % Initialize the segment diameter array (m)
    for seg = 1:Nseg
        %if (Ncell(seg) > 1)
        if (Ncell(seg) >= 1)
            D(seg) = Ncell(seg)*cell_size/pi;
        else
            D(seg) = 0;
        end
    end
    
    G = (pi*D.^4)./(128*mu*L);                                              % Calculate segment conductance
    
    for seg = 1:Nseg
        if (D(seg) ~= 0)
            H(seg) = (32*mu)/(pi*D(seg)^3);
        else
            H(seg) = 0;
        end
    end
    
    [P, Q, tau] = solve_for_flow_Ub1(G, Pin, Pout, H);                                % Solve for nodal pressure and segment flow
  
    plot_network_Ub1(segments, D, P, Q, tau, seg_cells)                            % Plot the network
    
    % Write to the video file
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame)

    % Wait for user to click on figure to continue
    if (t ~= Nt)
        %waitforbuttonpress
        pause(0.01)
    end
end

close(vidObj);
