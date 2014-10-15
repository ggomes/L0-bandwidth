function [  ] = bw_tests(  )

add_dependencies()

test_gaussian

testtest

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function []=test_gaussian()

A = load_SanPablo(100,0,'pretimed');
A.remove_intersection('Allston'); 
A.remove_intersection('Grayson');
A.remove_intersection('Dwight');
A.remove_intersection('Addison');
    
A.optimize();
A.plot();

disp(A.toString())

%%
function []=testtest()
% test
clear
close all

add_dependencies()

method = 'ip';
cycle = 70;
green_split = 0.5;
speed = 16;
numsamples = 500;
windowtype = 'pretimed';

A = class_artery(cycle,windowtype);
A.pretimed_optimization_method = method;
I = class_intersection('1',windowtype,green_split,green_split,0.185023324362067);
A.add_intersection(I);
A.add_segment(438,speed,speed);
I = class_intersection('2',windowtype,green_split,green_split,-0.026129390935360);
A.add_intersection(I);
A.initialize();

% optimize for pretimed
A.optimize();
A.plot();

pretimed.totalband = A.optbandwidth;
pretimed.bo = A.optbo;
pretimed.bi = A.optbi;
pretimed.reloffseto = [A.intersection.reloffseto];
pretimed.reloffset = [A.intersection.reloffseti];
pretimed.go = [A.intersection.go];
pretimed.gi = [A.intersection.gi];

% optimize for gaussian
A.setWindowType('gaussian')
A.optimize();
gaussian.totalband = A.optbandwidth;
gaussian.bo = A.optbo;
gaussian.bi = A.optbi;
gaussian.reloffseto = [A.intersection.reloffseto];
gaussian.reloffseti = [A.intersection.reloffseti];
gaussian.sigma_o = [A.intersection.sigma_o];
gaussian.sigma_i = [A.intersection.sigma_i];
gaussian.gamma_o = [A.intersection.gamma_o];
gaussian.gamma_i = [A.intersection.gamma_i];

% outbound, gaussian relative offset, pretimed green, pretimed bandwidth function
reloffset = gaussian.reloffseto;
green = pretimed.go;
mixed.bo = compute_bandwidth_pretimed(reloffset',green');

% inbound, gaussian relative offset, pretimed green, pretimed bandwidth function
reloffset = gaussian.reloffseti;
green = pretimed.gi;
mixed.bi = compute_bandwidth_pretimed(reloffset',green');

% plot
A.intersection(1).reloffseto = gaussian.reloffseto(1);
A.intersection(2).reloffseto = gaussian.reloffseto(2);
A.intersection(1).reloffseti = gaussian.reloffseti(1);
A.intersection(2).reloffseti = gaussian.reloffseti(2);
A.setWindowType('pretimed')
A.plot()

mixed.totalband = mixed.bo + mixed.bi;

pretimed.totalband-mixed.totalband

%%


