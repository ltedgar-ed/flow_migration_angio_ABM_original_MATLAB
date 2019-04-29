function [segments] = make_segments_Ub1(L)

v1 = [0; 1];
v2 = [1; 0];
v3 = [0; -1];
v4 = [0; 1];
v5 = [1; 0];
v6 = [0; -1];

vessel1 = zeros(6,2);

for seg = 1:5
    vessel1(seg+1,1) = 0;
    vessel1(seg+1,2) = sum(L(1:seg))*1e6;
end

vessel2 = zeros(11,2);
vessel2(1,:) = vessel1(6,:);

for seg = 6:15
    vessel2(seg-5+1,1) = sum(L(6:seg))*1e6;
    vessel2(seg-5+1,2) = 50;
end

vessel3 = zeros(5,2);
vessel3(1,:) = vessel2(11,:);

for seg = 16:20
    vessel3(seg-15+1,1) = 100;
    vessel3(seg-15+1,2) = 50 - sum(L(16:seg))*1e6;
end

vessel4 = zeros(6,2);
vessel4(1,:) = vessel1(6,:);

for seg = 21:25
    vessel4(seg-20+1,1) = 0;
    vessel4(seg-20+1,2) = 50 + sum(L(21:seg))*1e6;
end

vessel5 = zeros(11,2);
vessel5(1,:) = vessel4(6,:);

for seg = 26:35
    vessel5(seg-25+1,1) = sum(L(26:seg))*1e6;
    vessel5(seg-25+1,2) = 100;
end

vessel6 = zeros(5,2);
vessel6(1,:) = vessel5(11,:);

for seg = 36:40
    vessel6(seg-35+1,1) = 100;
    vessel6(seg-35+1,2) = 100 - sum(L(36:seg))*1e6;
end

segments1 = [vessel1; vessel2; vessel3];
segments2 = [vessel4; vessel5; vessel6];

segments = [segments1 segments2];