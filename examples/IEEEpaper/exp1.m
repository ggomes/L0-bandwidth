function [exp1name]=exp1(average_segment_length,cycle,green_split,max_abs_delta,speed,numsamples)
% Experiment #1
% -------------
% fix: cycle, green splits, delta splits
% independent variable: number of intersections
% random variable: Poisson segment lengths

add_dependencies()

rootfolder = fileparts(fileparts(mfilename('fullpath')));

exp1name = fullfile(rootfolder,'_storage',...
    ['exp1_'  ...
    num2str(average_segment_length) '_' ...
    num2str(cycle) '_' ...
    num2str(green_split*100,'%.2d') '_' ...
    num2str(max_abs_delta*100,'%.2d') '_' ...
    num2str(speed) '_' ...
    num2str(numsamples) ] );

if(exist([exp1name '.mat'],'file'))
    return
end

num_intersections = 2:8;
i=0;

% allocate
n_runs = length(num_intersections);
n_inter = num_intersections(end);

pretimed_simplex.totalband = nan(n_runs,numsamples);
pretimed_simplex.bo = nan(n_runs,numsamples);
pretimed_simplex.bi = nan(n_runs,numsamples);
pretimed_simplex.reloffseto = nan(n_inter,n_runs,numsamples);
pretimed_simplex.reloffseti = nan(n_inter,n_runs,numsamples);
pretimed_simplex.go = nan(n_inter,n_runs,numsamples);
pretimed_simplex.gi = nan(n_inter,n_runs,numsamples);

pretimed_ip.totalband = nan(n_runs,numsamples);
pretimed_ip.bo = nan(n_runs,numsamples);
pretimed_ip.bi = nan(n_runs,numsamples);
pretimed_ip.reloffseto = nan(n_inter,n_runs,numsamples);
pretimed_ip.reloffseti = nan(n_inter,n_runs,numsamples);
pretimed_ip.go = nan(n_inter,n_runs,numsamples);
pretimed_ip.gi = nan(n_inter,n_runs,numsamples);

gaussian.totalband = nan(n_runs,numsamples);
gaussian.bo = nan(n_runs,numsamples);
gaussian.bi = nan(n_runs,numsamples);
gaussian.reloffseto = nan(n_inter,n_runs,numsamples);
gaussian.reloffseti = nan(n_inter,n_runs,numsamples);
gaussian.sigma_o = nan(n_inter,n_runs,numsamples);
gaussian.sigma_i = nan(n_inter,n_runs,numsamples);
gaussian.gamma_o = nan(n_inter,n_runs,numsamples);
gaussian.gamma_i = nan(n_inter,n_runs,numsamples);

for n = num_intersections
    
    i=i+1;
    
    for j=1:numsamples
        
        % make the arterial
        windowtype = 'pretimed';
        A = class_artery(cycle,windowtype);
        for k=1:n-1
            delta = max_abs_delta*(2*rand(1)-1);
            I = class_intersection(num2str(k),windowtype,green_split,green_split,delta);
            A.add_intersection(I);
            A.add_segment(poissrnd(average_segment_length),speed,speed);
        end
        delta = max_abs_delta*(2*rand(1)-1);
        I = class_intersection(num2str(n),windowtype,green_split,green_split,delta);
        A.add_intersection(I);
        A.initialize();
        
        %             % optimize for pretimed simplex
        %             A.pretimed_optimization_method = 'simplex';
        %             A.optimize();
        %             pretimed_simplex.totalband(i,j) = A.optbandwidth;
        %             pretimed_simplex.bo(i,j) = A.optbo;
        %             pretimed_simplex.bi(i,j) = A.optbi;
        %             pretimed_simplex.reloffseto(1:n,i,j) = [A.intersection.reloffseto];
        %             pretimed_simplex.reloffseti(1:n,i,j) = [A.intersection.reloffseti];
        %             pretimed_simplex.go(1:n,i,j) = [A.intersection.go];
        %             pretimed_simplex.gi(1:n,i,j) = [A.intersection.gi];
        
        % optimize for pretimed ip
        A.pretimed_optimization_method = 'ip';
        A.optimize();
        pretimed_ip.totalband(i,j) = A.optbandwidth;
        pretimed_ip.bo(i,j) = A.optbo;
        pretimed_ip.bi(i,j) = A.optbi;
        pretimed_ip.reloffseto(1:n,i,j) = [A.intersection.reloffseto];
        pretimed_ip.reloffseti(1:n,i,j) = [A.intersection.reloffseti];
        pretimed_ip.go(1:n,i,j) = [A.intersection.go];
        pretimed_ip.gi(1:n,i,j) = [A.intersection.gi];
        
        % optimize for gaussian
        A.setWindowType('gaussian')
        A.optimize();
        gaussian.totalband(i,j) = A.optbandwidth;
        gaussian.bo(i,j) = A.optbo;
        gaussian.bi(i,j) = A.optbi;
        gaussian.reloffseto(1:n,i,j) = [A.intersection.reloffseto];
        gaussian.reloffseti(1:n,i,j) = [A.intersection.reloffseti];
        gaussian.sigma_o(1:n,i,j) = [A.intersection.sigma_o];
        gaussian.sigma_i(1:n,i,j) = [A.intersection.sigma_i];
        gaussian.gamma_o(1:n,i,j) = [A.intersection.gamma_o];
        gaussian.gamma_i(1:n,i,j) = [A.intersection.gamma_i];
        
    end
    
end

save(exp1name,'average_segment_length','green_split','num_intersections',...
    'speed','cycle','delta','max_abs_delta','numsamples',...
    'n_inter','gaussian','n_runs','pretimed_simplex','pretimed_ip')

