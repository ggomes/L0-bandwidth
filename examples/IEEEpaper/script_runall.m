clear
close all

root = fileparts(fileparts(mfilename('fullpath')));

facility = 'canal';
windowtype = 'gaussian';

ALL.average_segment_length = [400,700,1000];    % [ft]
ALL.cycle = [70,85,100,115];                    % seconds
ALL.green_split = [0.5,0.7];
ALL.max_abs_delta = [0,0.2];
ALL.speed = 16;
ALL.numsamples = 500;
save(fullfile(root,'_storage','ALL'),'ALL')

% experiment 0
cycle = 100;
delta = 0;
A=exp0(nan,facility,cycle,delta,'pretimed');

delta = 10/cycle;
A=exp0(nan,facility,cycle,delta,'pretimed');

Z = [...
[A.intersection.absoffseti]' ...
[A.intersection.absoffseto]' ... 
[A.intersection.reloffseti]' ... 
[A.intersection.reloffseto]' ]

% experiment 1
runnumber = 0;
for average_segment_length = ALL.average_segment_length
    for cycle = ALL.cycle
        for green_split = ALL.green_split
            for max_abs_delta = ALL.max_abs_delta
                for speed = ALL.speed
                    for numsamples = ALL.numsamples
                        runnumber = runnumber+1;
                        disp(['********** RUN NUMBER ' num2str(runnumber) ' ******************'])
                        exp1(average_segment_length,cycle,green_split,max_abs_delta,speed,numsamples);
                    end
                end
            end
        end
    end
end

% experiment 3
runnumber = 0;
for average_segment_length = ALL.average_segment_length
    for cycle = ALL.cycle
        for green_split = ALL.green_split
            for max_abs_delta = ALL.max_abs_delta
                for speed = ALL.speed
                    for numsamples = ALL.numsamples
                        runnumber = runnumber+1;
                        disp(['********** RUN NUMBER ' num2str(runnumber) ' ******************'])
                        exp3(average_segment_length,cycle,green_split,max_abs_delta,speed,numsamples);
                    end
                end
            end
        end
    end
end

% experiment 4
delta = 0;
exp4(facility,delta,windowtype)


