clear
close all

cycle = 100;
delta = 0.1;
windowtype = 'pretimed';

A = load_Canal(cycle,delta,windowtype);
A.pretimed_optimization_method = 'milp_ip';
tic
A.optimize();
toc
A.optbandwidth

[to,ti]=segment_travel_times(A);

A.modhalf(to)/A.cycle

B = load_Canal(cycle,delta,windowtype);
B.pretimed_optimization_method = 'ip';
tic
B.optimize();
toc
B.optbandwidth/B.cycle
  


