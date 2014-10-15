function []=exp2(exp1name)
% Experiment #2
% -------------
% Take the result of experiment #1,
% compute the bandwidith for pretimed windows resulting from the optimal
% offsets for gaussian windows.
% Compare that to the optimal pretimed bandwidths.

% Conclusion: the gaussian method also delivers to optimal for pretimed
% windows.

add_dependencies()

if(~exist([exp1name '.mat'],'file'))
    disp('run exp1 first')
    return
end

load(exp1name)

mixed_ip = compute_mixed(gaussian,pretimed_ip,n_runs,numsamples,num_intersections);
mixed_simplex = compute_mixed(gaussian,pretimed_simplex,n_runs,numsamples,num_intersections);

[ppt,op]=openppt('plots_exp_2',true);

figure('Visible','off')
subplot(211)
plot(pretimed_simplex.totalband'-mixed_simplex.totalband','.')
ylabel('simplex');
subplot(212)
plot(pretimed_ip.totalband'-mixed_ip.totalband','.')
ylabel('ip');

tit = ['seg_length=' num2str(average_segment_length) ' ' ...
    'cycle=' num2str(cycle) ' ' ...
    'green_split=' num2str(green_split*100,'%.2d') ' ' ...
    'delta_max=' num2str(max_abs_delta*100,'%.2d') ' ' ...
    'speed=' num2str(speed) ' ' ...
    'numsamples=' num2str(numsamples) ];

addslide(op,tit)
close

closeppt(ppt,op)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [mixed]=compute_mixed(gaussian,pretimed,n_runs,numsamples,num_intersections)

for i=1:n_runs
    for j=1:numsamples
        
        num_inter = num_intersections(i);
        
        % outbound, gaussian relative offset, pretimed green, pretimed bandwidth function
        reloffset = gaussian.reloffseto(1:num_inter,i,j);
        green = pretimed.go(1:num_inter,i,j);
        mixed.bo(i,j) = compute_bandwidth_pretimed(reloffset,green);
        
        % inbound, gaussian relative offset, pretimed green, pretimed bandwidth function
        reloffset = gaussian.reloffseti(1:num_inter,i,j);
        green = pretimed.gi(1:num_inter,i,j);
        mixed.bi(i,j) = compute_bandwidth_pretimed(reloffset,green);
    end
    
end

mixed.totalband = mixed.bo + mixed.bi;
