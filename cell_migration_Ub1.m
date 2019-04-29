% Subfunction to handle cellular migration in the ABM model
% subfunction of ABM_EC_polarize_migrate_flow_H.m
% 09/14/2018
% --------------------
function [seg_cells, new_seg_cells] = cell_migration_Ub1(seg, seg_cells, new_seg_cells, migrate, Q, tau, branch_rule, alpha)

% branch_rule = 2;

cell_size = 10*1e-6;                                                        % Set the size of each cell (m)

% Calculate chance of migrating depending on flow
%mchance = (1 - abs(Q(seg))/Qmax);                                   % Chance of a cell migrating depending on segment flow
mchance = 1;

if (mchance < 0)
    mchance = 0;
end

if (seg_cells{seg,1} ~= 0)
    for cell = 1:seg_cells{seg,1}                                   % Iterate through each cell of each segment and use a random number to determine if it migrates against flow
        mcell = rand();
        
        if (mcell <= mchance)                                       % If random number (0,1) is less than migration chacne...
            polar_vect = seg_cells{seg,2}(:,cell);
            migrate_vect = cell_size*polar_vect;
            
            % For cells in the first or fourth vertical vessel...
            if ((seg <= 5) || ((seg >= 21) && (seg <= 25)))                                          
                if (migrate_vect(2,1) >= cell_size/2)                   % If polarity vector is sufficiently aligned upwards, migrate upstream
                    seg_cells{seg,3}(1,cell) = 1;
                    migrate(seg) = migrate(seg) + 1;
                else if (migrate_vect(2,1) <= -cell_size/2)             % If polarity vector is sufficiently aligned downwards, migrate downstream
                        seg_cells{seg,3}(1,cell) = -1;
                        migrate(seg) = migrate(seg) + 1;
                    end
                end
            end
            
            % For cells in the second or fifth horizontal vessel...
            if (((seg >= 6) && (seg <= 15)) || ((seg >= 26) && (seg <= 35)))                         
                if (migrate_vect(1,1) >= cell_size/2)                   % If polarity vector is sufficiently aligned to the right, migrate upstream
                    seg_cells{seg,3}(1,cell) = 1;
                    migrate(seg) = migrate(seg) + 1;
                else if (migrate_vect(1,1) <= -cell_size/2)             % If polarity vector is sufficiently aligned to the left, migrate downstream
                        seg_cells{seg,3}(1,cell) = -1;
                        migrate(seg) = migrate(seg) + 1;
                    end
                end
            end
             
            % For cells in the third or sixth vertical vessel...
            if (((seg >= 16) && (seg <= 20)) || (seg >= 36))                         
                if (migrate_vect(2,1) >= cell_size/2)                   % If polarity vector is sufficiently aligned upwards, migrate downstream
                    seg_cells{seg,3}(1,cell) = -1;
                    migrate(seg) = migrate(seg) + 1;
                else if (migrate_vect(2,1) <= -cell_size/2)             % If polarity vector is sufficiently aligned downwards, migrate upstream
                        seg_cells{seg,3}(1,cell) = 1;
                        migrate(seg) = migrate(seg) + 1;
                    end
                end
            end
         
        end
    end
end

% Diffusion scheme to even out cell numbers in vertical vessel segments
down_seg = [];

if ((seg ~= 20) && (seg ~= 40))
    down_seg = seg + 1;
end

if (seg == 20)
    down_seg = 1;
end

if (seg == 40)
    if ((seg_cells{15,1} < seg_cells{16,1}) && (seg_cells{15,1} ~= 0))
        down_seg = 15;
    else
        down_seg = 16;
    end
end

if (seg == 5)
    if ((seg_cells{6,1} < seg_cells{21,1}) && (seg_cells{6,1} ~= 0))
        down_seg = 6;
    else
        down_seg = 21;
    end
end

if (isempty(down_seg) == 0)
    if (seg_cells{seg,1} ~= 0)
        if ((seg_cells{down_seg,1} < seg_cells{seg,1}) && (seg_cells{down_seg,1} ~= 0))
            seg_cells{seg,3}(randi([1 seg_cells{seg,1}])) = 0;
            migrate(seg) = migrate(seg) - 1;
        end
    end
end

% Move all cells that have been selected to migrate
if (migrate(seg) ~= 0)
    for cell = 1:seg_cells{seg,1}
        
        % If cell is migrating downstream (i.e., with flow if flow is positive)
        if (seg_cells{seg,3}(1,cell) == 1)
            new_seg_cells{seg,1} = new_seg_cells{seg,1} - 1;                % Decrease the number of cells in the parent segment
            cell_vect = seg_cells{seg,2}(:,cell);                           % Obtain the cell's polarity vector
            new_seg_cells{seg,2}(:,cell) = [0; 0];
            target = 0;
            
            % Move the cell and it's polarity vector to the target segment (downstream segment)
            if (seg ~= 20)
                target = seg+1;
            end
                       
            if (seg == 20)
                target = 1;
            end          
            
            if (seg == 40)
                target = 16;
            end
            
            switch branch_rule
                % Flow magnitude rule
                case 1
                    if (seg == 5)
                        if (Q(6) > Q(21))
                            target = 6;
                        else
                            target = 21;
                        end
                    end

                    if (seg == 40)
                        if (Q(16) > Q(15))
                            target = 16;
                        else
                            target = 15;
                        end
                    end
    
                % Polarity direction rule
                case 2
                    if (seg == 5)
                        if (dot(cell_vect, [0; 1]) > dot(cell_vect, [1; 0]))
                            target = 21;
                        else
                            target = 6;
                        end
                    end
                
                % Random rule
                case 3
                    if (seg == 5)
                        if (sign(rand() - 0.5) < 0)
                            target = 21;
                        else
                            target = 6;
                        end
                    end
                    
                case 4
                    if (seg == 5)
                        if (sign(rand() - 0.3) < 0)
                            target = 21;
                        else
                            target = 6;
                        end
                    end
                    
                case 5
                    Q1 = Q(6);
                    Q2 = Q(21);
                    Qratio = Q2/(Q1 + Q2);
                    
                    if (seg == 5)
                        if (sign(rand() - Qratio) < 0)
                            target = 21;
                        else
                            target = 6;
                        end
                    end
                    
                case 6
                    %alpha = 0.5;
                    
                    Q1 = Q(6) + Q(7) + Q(8);
                    Q2 = Q(21) + Q(22) + Q(23);
                    Qratio = Q2/(Q1 + Q2);
                    
%                     Nparent = (seg_cells{1,1} + seg_cells{2,1} + seg_cells{3,1} + seg_cells{4,1} + seg_cells{5,1});
%                     N2 = (seg_cells{21,1} + seg_cells{22,1} + seg_cells{23,1} + seg_cells{24,1} + seg_cells{25,1});
%                     N1 = (seg_cells{6,1} + seg_cells{7,1} + seg_cells{8,1} + seg_cells{9,1} + seg_cells{10,1});
%                     Nratio = N2/(N1 + N2);

                    Nparent = (seg_cells{3,1} + seg_cells{4,1} + seg_cells{5,1});
                    N2 = (seg_cells{21,1} + seg_cells{22,1} + seg_cells{23,1});
                    N1 = (seg_cells{6,1} + seg_cells{7,1} + seg_cells{8,1});
                    Nratio = N2/(N1 + N2);
                    
                    P2 = alpha*Qratio + (1 - alpha)*Nratio;
                    
                    if (seg == 5)
                        if (sign(rand() - P2) < 0)
                            target = 21;
                        else
                            target = 6;
                        end
                    end
                    
               case 7
                    tau1 = tau(6) + tau(7) + tau(8);
                    tau2 = tau(21) + tau(22) + tau(23);
                    tau_ratio = tau2/(tau1 + tau2);
                    
%                     Nparent = (seg_cells{1,1} + seg_cells{2,1} + seg_cells{3,1} + seg_cells{4,1} + seg_cells{5,1});
%                     N2 = (seg_cells{21,1} + seg_cells{22,1} + seg_cells{23,1} + seg_cells{24,1} + seg_cells{25,1});
%                     N1 = (seg_cells{6,1} + seg_cells{7,1} + seg_cells{8,1} + seg_cells{9,1} + seg_cells{10,1});
%                     Nratio = N2/(N1 + N2);

                    Nparent = (seg_cells{3,1} + seg_cells{4,1} + seg_cells{5,1});
                    N2 = (seg_cells{21,1} + seg_cells{22,1} + seg_cells{23,1});
                    N1 = (seg_cells{6,1} + seg_cells{7,1} + seg_cells{8,1});
                    Nratio = N2/(N1 + N2);
                    
                    P2 = alpha*tau_ratio + (1 - alpha)*Nratio;
                    
                    if (seg == 5)
                        if (sign(rand() - P2) < 0)
                            target = 21;
                        else
                            target = 6;
                        end
                    end
            end
            
                
            
            % Rotate polarity vectors of cells crossing into a branch
            if ((seg == 5) && (target == 6))
                theta = -90;
                cell_vect = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]*cell_vect;
                cell_vect = cell_vect/norm(cell_vect);
            end
            
            if ((seg == 25) && (target == 26))
                theta = -90;
                cell_vect = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]*cell_vect;
                cell_vect = cell_vect/norm(cell_vect);
            end
            
            if ((seg == 15) && ((target == 16) || (target == 40)))
                theta = -90;
                cell_vect = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]*cell_vect;
                cell_vect = cell_vect/norm(cell_vect);
            end
              
            if ((seg == 35) && (target == 36))
                theta = -90;
                cell_vect = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]*cell_vect;
                cell_vect = cell_vect/norm(cell_vect);
            end
            
            if ((seg == 20) && (target == 1))
                theta = -180;
                cell_vect = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]*cell_vect;
                cell_vect = cell_vect/norm(cell_vect);
            end
            
            new_seg_cells{target,1} = new_seg_cells{target,1} + 1;
            new_seg_cells{target,2} = [new_seg_cells{target,2} cell_vect];            
        end
        
        % If cell is migrating upstream (i.e., against flow if flow is positive)...
        if (seg_cells{seg,3}(1,cell) == -1)
            new_seg_cells{seg,1} = new_seg_cells{seg,1} - 1;                % Decrease the number of cells in the parent segment
            cell_vect = seg_cells{seg,2}(:,cell);                           % Obtain the cell's polarity vector
            new_seg_cells{seg,2}(:,cell) = [0; 0];
            target = 0;
            
            % Move the cell and it's polarity vector to the target segment (upstream segment)
            if (seg ~= 1)
                target = seg-1;
            end
            
            if (seg == 1)
                target = 20;
            end          
            
            if (seg == 21)
                target = 5;
            end
            
            switch branch_rule
                % Flow magnitude rule
                case 1
                    if (seg == 16)
                        if (Q(15) > Q(40))
                            target = 15;
                        else
                            target = 40;
                        end
                    end

                    if (seg == 6)
                        if (Q(5) > Q(21))
                            target = 5;
                        else
                            target = 21;
                        end
                    end

                % Polarity direction rule
                case 2
                    if (seg == 16)
                        if (dot(cell_vect, [0; 1]) > dot(cell_vect, [-1; 0]))
                            target = 40;
                        else
                            target = 15;
                        end
                    end
                
                % Random rule
                case 3
                    if (seg == 16)
                        if (sign(rand() - 0.5) < 0)
                            target = 40;
                        else
                            target = 15;
                        end
                    end

                case 4
                    if (seg == 16)
                        if (sign(rand() - 0.3) < 0)
                            target = 40;
                        else
                            target = 15;
                        end
                    end
                case 5
                    Qratio = Q(40)/(Q(40) + Q(15));
                    
                    if (seg == 16)
                        if (sign(rand() - Qratio) < 0)
                            target = 40;
                        else
                            target = 15;
                        end
                    end
                    
                case 6
                    %alpha = 0.5;
                    
                    Q1 = Q(15) + Q(14) + Q(13);
                    Q2 = Q(40) + Q(39) + Q(38);
                    Qratio = Q2/(Q1 + Q2);
                    
%                     Nparent = (seg_cells{16,1} + seg_cells{17,1} + seg_cells{18,1} + seg_cells{19,1} + seg_cells{20,1});
%                     N2 = (seg_cells{40,1} + seg_cells{39,1} + seg_cells{38,1} + seg_cells{37,1} + seg_cells{36,1});
%                     N1 = (seg_cells{15,1} + seg_cells{14,1} + seg_cells{13,1} + seg_cells{12,1} + seg_cells{11,1});
%                     Nratio = N2/(N1 + N2);

                    Nparent = (seg_cells{16,1} + seg_cells{17,1} + seg_cells{18,1});
                    N2 = (seg_cells{40,1} + seg_cells{39,1} + seg_cells{38,1});
                    N1 = (seg_cells{15,1} + seg_cells{14,1} + seg_cells{13,1});
                    Nratio = N2/(N1 + N2);
                    
                    P2 = alpha*Qratio + (1 - alpha)*Nratio;
                    
                    if (seg == 16)
                        if (sign(rand() - P2) < 0)
                            target = 40;
                        else
                            target = 15;
                        end
                    end
                    
               case 7
                    tau1 = tau(15) + tau(14) + tau(13);
                    tau2 = tau(40) + tau(39) + tau(38);
                    tau_ratio = tau2/(tau1 + tau2);
                    
%                     Nparent = (seg_cells{16,1} + seg_cells{17,1} + seg_cells{18,1} + seg_cells{19,1} + seg_cells{20,1});
%                     N2 = (seg_cells{40,1} + seg_cells{39,1} + seg_cells{38,1} + seg_cells{37,1} + seg_cells{36,1});
%                     N1 = (seg_cells{15,1} + seg_cells{14,1} + seg_cells{13,1} + seg_cells{12,1} + seg_cells{11,1});
%                     Nratio = N2/(N1 + N2);
                    
                    Nparent = (seg_cells{16,1} + seg_cells{17,1} + seg_cells{18,1});
                    N2 = (seg_cells{40,1} + seg_cells{39,1} + seg_cells{38,1});
                    N1 = (seg_cells{15,1} + seg_cells{14,1} + seg_cells{13,1});
                    Nratio = N2/(N1 + N2);
                    
                    P2 = alpha*tau_ratio + (1 - alpha)*Nratio;
                    
                    if (seg == 16)
                        if (sign(rand() - P2) < 0)
                            target = 40;
                        else
                            target = 15;
                        end
                    end
            end           

            % Rotate polarity vectors of cells crossing into a branch  
            if (((seg == 16) || (seg == 40)) && (target == 15))
                theta = 90;
                cell_vect = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]*cell_vect;
                cell_vect = cell_vect/norm(cell_vect);
            end
            
            if ((seg == 36) && (target == 35))
                theta = 90;
                cell_vect = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]*cell_vect;
                cell_vect = cell_vect/norm(cell_vect);
            end
            
            if ((seg == 26) && (target == 25))
                theta = 90;
                cell_vect = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]*cell_vect;
                cell_vect = cell_vect/norm(cell_vect);
            end
            
            if ((seg == 6) && ((target == 5) || (target == 21)))
                theta = 90;
                cell_vect = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]*cell_vect;
                cell_vect = cell_vect/norm(cell_vect);
            end
              
            if ((seg == 1) && (target == 20))
                theta = 180;
                cell_vect = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]*cell_vect;
                cell_vect = cell_vect/norm(cell_vect);
            end
            
            new_seg_cells{target,1} = new_seg_cells{target,1} + 1;
            new_seg_cells{target,2} = [new_seg_cells{target,2} cell_vect];
        end
    end
end
        