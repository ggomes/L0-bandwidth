function []=exp3(average_segment_length,cycle,green_split,max_abs_delta,speed,numsamples)
% Experiment #3
% -------------
% Take the result of experiment #1,
% For each of the 7x500 runs, do a series of perturbations; add and
% subtract 1 second from each of the offsets and compute resulting pretimed
% bandwidths. Do this for pretimed and gaussian offsets separately. Then
% compare the perturbed bandwidths to show that the gaussian solution is
% more robust.
% But in fact the pretimed seems more robust.

rootfolder = fileparts(fileparts(mfilename('fullpath')));

exp1name = fullfile(rootfolder,'_storage',...
    ['exp1_'  ...
    num2str(average_segment_length) '_' ...
    num2str(cycle) '_' ...
    num2str(green_split*100,'%.2d') '_' ...
    num2str(max_abs_delta*100,'%.2d') '_' ...
    num2str(speed) '_' ...
    num2str(numsamples)]);

if(~exist([exp1name '.mat'],'file'))
    disp('run exp1 first')
    return
end

exp3name = regexprep(exp1name,'exp1','exp3');
if(exist([exp3name '.mat'],'file'))
    return
end

load(exp1name)

pert_pretimed_simplex   = nan(n_runs,numsamples);
pert_pretimed_ip        = nan(n_runs,numsamples);
pert_gaussian           = nan(n_runs,numsamples);

for i=1:n_runs
    for j=1:numsamples
        
        num_inter = num_intersections(i);
        
        green_i = pretimed_ip.gi(1:num_inter,i,j);
        green_o = pretimed_ip.go(1:num_inter,i,j);
        totalband = pretimed_ip.totalband(i,j);
        
        % pretimed simplex ........................................................
        pert_pretimed_simplex(i,j) = compute_pert( ...
            pretimed_simplex.reloffseti(1:num_inter,i,j) , ...
            pretimed_simplex.reloffseto(1:num_inter,i,j) , ...
            totalband,num_inter,green_o,green_i,cycle);

        % pretimed simplex ........................................................
        pert_pretimed_ip(i,j) = compute_pert( ...
            pretimed_ip.reloffseti(1:num_inter,i,j) , ...
            pretimed_ip.reloffseto(1:num_inter,i,j) , ...
            totalband,num_inter,green_o,green_i,cycle);
        
        % gaussian ........................................................
        pert_gaussian(i,j) = compute_pert( ...
            gaussian.reloffseti(1:num_inter,i,j) , ...
            gaussian.reloffseto(1:num_inter,i,j) , ...
            totalband,num_inter,green_o,green_i,cycle);
    end
    
end

save(exp3name,'exp1name','pert_pretimed_simplex','pert_pretimed_ip','pert_gaussian')

% -=-------------------------------------------------------------------
function [var]=compute_pert(reloffset_o,reloffset_i,band_nopert,num_inter,green_o,green_i,cycle)

pertsize = 1;
band_pert = [];
for k=1:num_inter
    
    % perturb +pertsize second
    reloffset_pert_o = reloffset_o;
    reloffset_pert_i = reloffset_i;
    reloffset_pert_o(k) = modhalf(reloffset_pert_o(k)+pertsize,cycle);
    reloffset_pert_i(k) = modhalf(reloffset_pert_i(k)+pertsize,cycle);
    
    % compute pertubed bandwidth
    band_pert_o = compute_bandwidth_pretimed(reloffset_pert_o,green_o);
    band_pert_i = compute_bandwidth_pretimed(reloffset_pert_i,green_i);
    band_pert = [band_pert band_pert_o+band_pert_i];
    
    % perturb -pertsize second
    reloffset_pert_o = reloffset_o;
    reloffset_pert_i = reloffset_i;
    reloffset_pert_o(k) = modhalf(reloffset_pert_o(k)-pertsize,cycle);
    reloffset_pert_i(k) = modhalf(reloffset_pert_i(k)-pertsize,cycle);
    
    % compute pertubed bandwidth
    band_pert_o = compute_bandwidth_pretimed(reloffset_pert_o,green_o);
    band_pert_i = compute_bandwidth_pretimed(reloffset_pert_i,green_i);
    band_pert = [band_pert band_pert_o+band_pert_i];
    
end

var = sqrt(mean((band_pert - band_nopert).^2));
