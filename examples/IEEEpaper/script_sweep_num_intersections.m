clear
close all

average_segment_length = 1500;      % [ft]
cycle = 60;                        % seconds
green_split = 0.5;
v_in = 30;
v_out = 30;
max_abs_delta = 0.25;

num_intersections = 2:8;
numsamples = 10;
i=0;
for n = num_intersections
    
    i=i+1;

    for j=1:numsamples

        [i j]

        % make the arterial
        windowtype = 'pretimed';
        A = class_artery(cycle,windowtype);

        for k=1:n-1
            delta = max_abs_delta*(2*rand(1)-1);
            I = class_intersection(num2str(k),windowtype,green_split,green_split,delta);
            A.add_intersection(I);
            A.add_segment(poissrnd(average_segment_length),v_in,v_out);
        end
        delta = max_abs_delta*(2*rand(1)-1);
        I = class_intersection(num2str(num_intersections),windowtype,green_split,green_split,delta);
        A.add_intersection(I);
        A.initialize();

        % optimize for pretimed
        A.optimize();
        B_pretimed(i,j) = A.compute_total_bandwidth();

        AllAs{i,j} = A;
        clear A
        %A.plot();
        
        % optimize for gaussian
%         A.setWindowType('gaussian')
%         A.optimize();
%         B_gaussian(i,j) = A.compute_total_bandwidth();

    end
    
end

figure
plot(B_pretimed)

% figure
% plot(B_gaussian)

for i=1:size(AllAs,1)
    for j=1:size(AllAs,2)
        optbandwidth(i,j) = AllAs{1,i}.optbandwidth;
    end
end
