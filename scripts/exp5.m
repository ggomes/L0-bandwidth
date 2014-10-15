function []=exp5(ALL)
% Experiment #5
% -------------
% Surface plots

add_dependencies()

speed = ALL.speed(1);
numsamples = ALL.numsamples(1);

n1 = length(ALL.average_segment_length);
n2 = length(ALL.cycle);
n3 = length(ALL.green_split);
n4 = length(ALL.max_abs_delta);

PretimedBand = nan(n1,n2,n3,n4,7);
GaussianBand = nan(n1,n2,n3,n4,7);

for i1 = 1:n1
    for i2 = 1:n2
        for i3 = 1:n3
            for i4 = 1:n4
                
                average_segment_length = ALL.average_segment_length(i1);
                cycle = ALL.cycle(i2);
                green_split = ALL.green_split(i3);
                max_abs_delta = ALL.max_abs_delta(i4);
                
                exp1name = ['exp1_'  ...
                    num2str(average_segment_length) '_' ...
                    num2str(cycle) '_' ...
                    num2str(green_split*100,'%.2d') '_' ...
                    num2str(max_abs_delta*100,'%.2d') '_' ...
                    num2str(speed) '_' ...
                    num2str(numsamples) ];
                
                
                load(exp1name)
                
                PretimedBand(i1,i2,i3,i4,:) = mean(pretimed_ip.totalband,2);
                GaussianBand(i1,i2,i3,i4,:) = mean(gaussian.totalband,2);
            end
        end
    end
end
clear pretimed_* gaussian


figure
for i1 = 1:n1
    for i3 = 1:n3
        for i4 = 1:n4
            
            average_segment_length = ALL.average_segment_length(i1); 
            %cycle = ALL.cycle(i2);
            green_split = ALL.green_split(i3);
            max_abs_delta = ALL.max_abs_delta(i4);
      
            clf
            plot( reshape(PretimedBand(i1,:,i3,i4,:),n2,7)' )
            legend(num2str(ALL.cycle'))
            title(['length=' num2str(average_segment_length) ' green_split=' num2str(green_split) ' delta=' num2str(max_abs_delta)])
            
            pause
        end
        
    end
end





