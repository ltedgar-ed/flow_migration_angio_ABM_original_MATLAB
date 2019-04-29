% Solve for the flow in the bifuracting vessel network
% subfunction of simple_ABM_EC_migration_flow.m
% 08/08/2018
% --------------------
function [P, Q, tau] = solve_for_flow_Ub1(G, Pin, Pout, H)

P1 = Pin;                                          % Inlet pressure 1 (Pa)
P21 = Pout;                                        % Outlet pressure 2 (Pa)
    
Nn = 40;                                            % Number of nodes
Nseg = 40;                                          % Number of segments

% Set any zero values for conductance to really small value
for seg = 1:Nseg
    if (G(seg) == 0)
        G(seg) = 1e-25;
    end
end

P = zeros(Nn,1);                                    % Initialize the nodal pressure array (Pa)
Q = zeros(Nseg,1);                                  % Initialize the segment flow array (m3/s)

C = zeros(Nn,Nn);                                   % Initialize the conductance matrix (m4/Pa-s-m)
B = zeros(Nn,1);                                    % Intialize the solution array (m3/s)

% Set equation for node 1
C(1,1) = G(1)*1;
B(1,1) = G(1)*P1;

% Set equation for node 2
C(2,:) = [0 (G(1)+G(2)) -G(2) zeros(1,Nn-3)];
B(2,1) = G(1)*P1;

% Set equations for nodes 3 through 5
for seg = 3:5
    C(seg,:) = [zeros(1,seg-2) -G(seg-1) (G(seg-1)+G(seg)) -G(seg) zeros(1,Nn-seg-1)];
end

% Set equation for node 6
%C(6,:) = [zeros(1,4) -G(5) (G(5) + G(6) + G(21)) -G(6) zeros(1,14) -G(21) zeros(1,18)];
C(6,5) = -G(5);
C(6,6) = (G(5) + G(6) + G(21));
C(6,7) = -G(6);
C(6,22) = -G(21);

% Set equations for nodes 7 through 15
for seg = 7:15
    C(seg,:) = [zeros(1,seg-2) -G(seg-1) (G(seg-1)+G(seg)) -G(seg) zeros(1,Nn-seg-1)];
end

% Set equation for node 16
C(16,15) = -G(15);
C(16,16) = (G(15) + G(16) + G(40));
C(16,17) = -G(16);
C(16,40) = -G(40);

% Set equations for nodes 17 through 19
for seg = 17:19
    C(seg,:) = [zeros(1,seg-2) -G(seg-1) (G(seg-1)+G(seg)) -G(seg) zeros(1,Nn-seg-1)];
end

% Set equation for node 20
%C(30,:) = [zeros(1,28) -G(29) (G(29)+G(30)) 0];
C(20,19) = -G(19);
C(20,20) = (G(19) + G(20));
B(20,1) = G(20)*P21;

% Set equation for node 21
C(21,21) = G(20)*1;
B(21,1) = G(20)*P21;

% Set equation for node 22
C(22,6) = -G(21);
C(22,22) = (G(21) + G(22));
C(22,23) = -G(22);

% Set equations for nodes 23 through 39
for seg = 23:39
    C(seg,:) = [zeros(1,seg-2) -G(seg-1) (G(seg-1)+G(seg)) -G(seg) zeros(1,Nn-seg-1)];
end

% Set equation for node 40
C(40,16) = -G(40);
C(40,39) = -G(39);
C(40,40) = (G(39) + G(40));

P = C\B;

% Calculate flow
for seg = 1:20
    Q(seg) = -G(seg)*(P(seg+1)-P(seg));
end

Q(21) = -G(21)*(P(22) - P(6));

for seg = 22:39
    Q(seg) = -G(seg)*(P(seg+1)-P(seg));
end

Q(40) = -G(40)*(P(16) - P(40));

for i = 1:length(G)
    tau(i) = H(i)*Q(i);
end
