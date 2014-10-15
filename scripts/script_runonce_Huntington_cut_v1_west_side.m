clear
close all

addpath(fullfile(fileparts(mfilename('fullpath')),'src'))

% decide of the cycle
% for now, east part, let's fix it at: (112+112+156+136)/4 = 516/4 = 129
% we could also decide of weights
cycle = 129;

windowtype = 'pretimed';

% let's fix delta at 0
% this has to be checked as delta = 20 on Huntington Dr & Santa Clara St
delta = 0;

A = load_Huntington_cut_v1_west_side(cycle,delta,windowtype);

clf

A.optimize();
A.plot(gcf);

disp('done')