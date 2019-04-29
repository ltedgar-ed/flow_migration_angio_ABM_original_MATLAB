% Realign the polarity vectors of all cells in a segment based on weights
% subfunction of ABM_EC_polarize_migrate_flow_H.m
% 09/14/2018
% --------------------
function [seg_cells, new_seg_cells] = realign_polarity_Ub1(seg, Q, seg_cells, new_seg_cells, w1, w2, w3, w4)

% Iterate through each cell in the segment and realign its polarity vector
if (seg_cells{seg,1} ~= 0)
    for cell = 1:seg_cells{seg,1}
        
        % Persistance alignment component
        polar_vect = seg_cells{seg,2}(:,cell);
        
        % Flow alignment component
        if (((seg >= 1) && (seg <= 5)) || ((seg >=21) && (seg <= 25)))
            flow_vect = sign(Q(seg))*[0; 1];
        end
        
        if (((seg >= 6) && (seg <= 15)) || ((seg >= 26) && (seg <= 35)))
            flow_vect = sign(Q(seg))*[1; 0];
        end
        
        if (((seg >= 16) && (seg <= 20)) || ((seg >= 36) && (seg <= 40))) 
            flow_vect = sign(Q(seg))*[0; -1];
        end
        
        flow_vect = -flow_vect;
        
%         % Alignment from neighbors
%         neighbors = [];
%         
%         if (seg == 1)
%             num_neigh = 2;
%             neighbors = [2 10];
%         end
%         
%         if (seg == 5)
%             num_neigh = 3;
%             neighbors = [4 6 11];
%         end
%         
%         if (seg == 6)
%             num_neigh = 3;
%             neighbors = [5 11 7];
%         end
%         
%         if (seg == 10)
%             num_neigh = 2;
%             neighbors = [1 9];
%         end
%         
%         if (seg == 11)
%             num_neigh = 3;
%             neighbors = [5 6 12];
%         end
%         
%         if (seg == 20)
%             num_neigh = 3;
%             neighbors = [19 25 26];
%         end
%         
%         if (seg == 21)
%             num_neigh = 2;
%             neighbors = [22 30];
%         end
%         
%         if (seg == 25)
%             num_neigh = 3;
%             neighbors = [20 24 26];
%         end
%         
%         if (seg == 26)
%             num_neigh = 3;
%             neighbors = [20 25 26];
%         end
%         
%         if (seg == 30)
%             num_neigh = 2;
%             neighbors = [29 21];
%         end
%         
%         if ((seg ~= 1) && (seg ~= 5) && (seg ~= 6) && (seg ~= 10) && (seg ~= 11) && (seg ~= 20) && (seg ~= 21) && (seg ~= 25) && (seg ~= 26) && (seg ~= 30))
%             num_neigh = 2;
%             neighbors = [seg-1 seg+1];
%         end
%         
%         cell_count = 0;
%         neigh_vect = [0; 0];
%         
%         for nseg = 1:num_neigh
%             neigh_seg = neighbors(nseg);
%             
%             if (seg_cells{neigh_seg,1} ~= 0)
%                 for neigh_cell = 1:seg_cells{neigh_seg,1}
%                     neigh_vect = neigh_vect + seg_cells{neigh_seg,2}(:,neigh_cell);
%                     cell_count = cell_count + 1;
%                 end
%             end
%         end
%         
%         neigh_vect = neigh_vect/cell_count;
%         neigh_vect = neigh_vect/norm(neigh_vect);
        
        % Calcualte the flow alignment angle
        dotprod2 = dot(flow_vect, polar_vect);
        
        if (abs(dotprod2) > 1)
            dotprod2 = sign(dotprod2)*1;
        end
        
        phi2 = acosd(dotprod2);
        
        
        cr2 = cross([polar_vect; 0], [flow_vect; 0]);
        
        if (cr2(3) < 0)
            phi2 = -phi2;
        end
               
%         % Calculate the neighboring polarity angle
%         dotprod3 = dot(neigh_vect, polar_vect);
%         
%         if (abs(dotprod3) > 1)
%             dotprod3 = sign(dotprod3)*1;
%         end
%         
%         phi3 = acosd(dotprod3);
%         
%         cr3 = cross([polar_vect; 0], [neigh_vect; 0]);
%         
%         if (cr3(3) < 0)
%             phi3 = -phi3;
%         end
%         
%         if (isnan(phi3))
%             phi3 = 0;
%         end
        
        % Calculate the random walk angle
        rand_walk_vect = randn(2,1);
        rand_walk_vect = rand_walk_vect/norm(rand_walk_vect);
        
        dotprod4 = dot(rand_walk_vect, polar_vect);
        
        if (abs(dotprod4) > 1)
            dotprod4 = sign(dotprod4)*1;
        end
        
        phi4 = acosd(dotprod4);
        
        cr4 = cross([polar_vect; 0], [rand_walk_vect; 0]);
        
        if (cr4(3) < 0)
            phi4 = -phi4;
        end
                
        % Calculate the new polarity vector
        %theta = w1*0 + w2*phi2 + w3*phi3 + w4*phi4;
        theta = w1*0 + w2*phi2 + w4*phi4;
        
        new_polar_vect = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]*polar_vect;
        new_polar_vect = new_polar_vect/norm(new_polar_vect);
        
        if ((isnan(new_polar_vect(1))) || (isnan(new_polar_vect(2))))
            disp('theta is NaN')
            pause
        end
        
        seg_cells{seg,2}(:,cell) = new_polar_vect;
        new_seg_cells{seg,2}(:,cell) = new_polar_vect;
    end
end