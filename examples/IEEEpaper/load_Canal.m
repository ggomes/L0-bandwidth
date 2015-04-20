function [A]=load_Canal(cycle,delta,windowtype)

A = class_artery(cycle,windowtype);

speed = 16;

g = 1-21/70;
I = class_intersection('1',windowtype,g,g,delta);
A.add_intersection(I);
A.add_segment(371,speed,speed);

g = 1-24/70;
I = class_intersection('2',windowtype,g,g,delta);
A.add_intersection(I);
A.add_segment(371,speed,speed);

g = 1-35/70;
I = class_intersection('3',windowtype,g,g,delta);
A.add_intersection(I);
A.add_segment(390,speed,speed);

g = 1-34/70;
I = class_intersection('4',windowtype,g,g,delta);
A.add_intersection(I);
A.add_segment(276,speed,speed);

g = 1-28/70;
I = class_intersection('5',windowtype,g,g,delta);
A.add_intersection(I);
A.add_segment(472,speed,speed);

g = 1-23/70;
I = class_intersection('6',windowtype,g,g,delta);
A.add_intersection(I);
A.add_segment(364,speed,speed);

g = 1-42/70;
I = class_intersection('7',windowtype,g,g,delta);
A.add_intersection(I);
A.add_segment(968,speed,speed);

g = 1-42/70;
I = class_intersection('8',windowtype,g,g,delta);
A.add_intersection(I);
A.add_segment(377,speed,speed);

g = 1-21/70;
I = class_intersection('9',windowtype,g,g,delta);
A.add_intersection(I);

A.initialize();