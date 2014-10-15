function [A] = load_Huntington_cut_v1_west_side(cycle,delta,windowtype)

% this function built a specific artery object

A = class_artery(cycle,windowtype);

% we decide to cut Huntington in two parts
% because of the different cycles
% we cut according to similar cycles in two parts
% cf load_Huntington_cut_v1_east_side

% arguments of classe_intersection:
% name
% windowtype?
% portion of green outbound
% portion of green inbound
% delta

% arguments of add_segment:
% length in feet
% speed in miles/hour in
% speed in miles/hour out

% speed limits gathered from http://www.google.com/mapmaker

% % cycle 66
% % delta 0
% % WB & EB speed limit bewteen 210 WB Off-Ramp & 210 EB Off-Ramp: 35
% % length bewteen 210 WB Off-Ramp & 210 EB Off-Ramp: 722
% I = class_intersection('210WBOR',windowtype,30/66,30/66,delta);
% A.add_intersection(I);
% A.add_segment(722,35,35);
% 
% % cycle 84
% % delta 0
% % WB & EB speed limit bewteen 210 EB Off-Ramp & 5th Ave: 35
% % length bewteen 210 EB Off-Ramp & 5th Ave: 437
% I = class_intersection('210EBOR',windowtype,36/84,36/84,delta);
% A.add_intersection(I);
% A.add_segment(437,35,35);
% 
% % cycle 88
% % delta 0
% % WB speed limit bewteen 5th Ave & Gateway Dr: 35
% % EB speed limit bewteen 5th Ave & Gateway Dr: 30
% % length bewteen 5th Ave & Gateway Dr: 597
% I = class_intersection('5th',windowtype,50/88,50/88,delta);
% A.add_intersection(I);
% A.add_segment(597,30,35);
% 
% % cycle 88
% % delta 0
% % WB & EB speed limit bewteen Gateway Dr & 2nd Ave: 30
% % length bewteen Gateway Dr & 2nd Ave: 560
% I = class_intersection('Gateway',windowtype,50/88,50/88,delta);
% A.add_intersection(I);
% A.add_segment(560,30,30);

% cycle 112
% delta 0
% WB & EB speed limit bewteen 2nd Ave & 1st Ave: 30
% length bewteen 2nd Ave & 1st Ave: 956
I = class_intersection('2nd',windowtype,50/112,50/112,delta);
A.add_intersection(I);
A.add_segment(956,30,30);

% cycle 112
% delta 0
% WB & EB speed limit bewteen 1st Ave & Santa Anita Ave: 30
% length bewteen 1st Ave & Santa Anita Ave: 837
I = class_intersection('1st',windowtype,50/112,50/112,delta);
A.add_intersection(I);
A.add_segment(837,30,30);

% cycle 156
% delta 0
% WB & EB speed limit bewteen Santa Anita Ave & Santa Clara St: 30
% length bewteen Santa Anita Ave & Santa Clara St: 998
I = class_intersection('SantaAnita',windowtype,50/156,50/156,delta);
A.add_intersection(I);
A.add_segment(998,30,30);

% cycle 136
% delta 20
I = class_intersection('SantaClara',windowtype,50/136,50/136,delta);
A.add_intersection(I);

A.initialize();