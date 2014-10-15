clear
close all

addpath(fullfile(fileparts(mfilename('fullpath')),'src'))

% decide of the cycle
% for now, east part, let's fix it at: (66+84+88+88)/4 = 326/4 = 81.5
% we could also decide of weights
cycle = 81.5;

windowtype = 'pretimed';

% let's fix delta at 0
delta = 0;

A = load_Huntington_cut_v1_east_side(cycle,delta,windowtype);

clf

A.optimize();
A.plot(gcf);

disp('done')