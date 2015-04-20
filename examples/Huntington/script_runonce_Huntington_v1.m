clear
close all

% decide of the cycle
% maybe we'll cut Huntington in two parts because of the different cycles
% for now, let's fix it at: (66+84+88+88+112+112+156+136)/8 = 842/8 = 105.25
% THIS IS AN ARBITRARY CHOICE, we could have decided of weights
cycle = 105;

windowtype = 'pretimed';

% see about delta, it depends on the part of Huntington
% for now, let's fix it at 0 so far as it is only nonnegative on Huntington Dr & Santa Clara St
delta = 0;

A = load_Huntington_v1(cycle,delta,windowtype);

clf

A.optimize();
A.plot(gcf);

disp('done')