function []=script_run_random()

numsamples = 100;

param.seg_length        = struct('mean',nan,'dev',nan,'min',100,'max',1100,'dom','N');
param.cycle             = struct('mean',nan,'dev',nan,'min',50,'max',180,'dom','N');
param.green_split       = struct('mean',0.5,'dev',0.3,'min',0.1,'max',0.7,'dom','R');
param.delta             = struct('mean',0.1,'dev',0.1,'min',0,'max',0.2,'dom','R');
param.speed             = struct('mean',nan,'dev',nan,'min',5,'max',25,'dom','N');
param.num_intersections = struct('mean',nan,'dev',nan,'min',1,'max',8,'dom','N');

add_dependencies()

% rootfolder = fileparts(fileparts(mfilename('fullpath')));

% allocate
X = repmat( struct('lp',[],'milp',[],'param',[],'lp_time',nan,'milp_time',nan) , 1,numsamples);
clear x

for j=1:numsamples
    
    % make the arterial
    windowtype = 'pretimed';
    P = sample(param);
    A = class_artery(P.cycle,windowtype);
    n = P.num_intersections;
    for k=1:n-1
        I = class_intersection(num2str(k),windowtype,P.green_split(k),P.green_split(k),P.delta(k));
        A.add_intersection(I);
        A.add_segment(P.seg_length(k),P.speed(k),P.speed(k));
    end
    I = class_intersection(num2str(n),windowtype,P.green_split(n),P.green_split(n),P.delta(n));
    A.add_intersection(I);
    A.initialize();
     
    X(j).param = P;
        
    % lp
    A.pretimed_optimization_method = 'ip';
    tic, A.optimize();
    X(j).lp_time = toc;
    X(j).lp = A;

    % milp
    A.pretimed_optimization_method = 'milp_ip';
    tic, A.optimize();
    X(j).milp_time = toc;
    X(j).milp = A;
       
    clear P A
end

save trial

Z = [X.lp];
lp.time = [X.lp_time];
lp.bw=[Z.optbandwidth];

Z = [X.milp];
milp.time = [X.milp_time];
milp.bw=[Z.optbandwidth];


figure('Position',[240   190   807   380])
plot(lp.bw,'ko','LineWidth',2,'MarkerSize',10)
hold on
plot(milp.bw,'r+','LineWidth',2,'MarkerSize',10)
grid
ylabel('optimal bandwidth')
xlabel('trial number')
legend('lp','milp')

figure
plot(lp.time,'k-o')
hold on
plot(milp.time,'r--+')
ylabel('execution time')
xlabel('trial number')


function [P]=sample(param)
F = fieldnames(param);
n = sample_value(param.num_intersections,1);
P.num_intersections = n;
for i = 1:length(F)
    name = F{i};
    switch name
        case {'cycle'}      % 1 value
            P.(F{i}) = sample_value(param.(F{i}),1);
        case {'seg_length','speed'} % n-1 value
            P.(F{i}) = sample_value(param.(F{i}),n-1);
        case {'delta','green_split'}     % n value
            P.(F{i}) = sample_value(param.(F{i}),n);
    end
end

function [v]=sample_value(z,n)
switch z.dom
    case 'R'
        v = z.mean + z.dev*2*(rand(1,n)-0.5);
        v = min([v;z.max*ones(1,n)]);
        v = max([v;z.min*ones(1,n)]);
    case 'N'
        v = z.min + randi(z.max-z.min+1,1,n) -1;
end
