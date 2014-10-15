function []=script_makeallplots()
add_dependencies()

load(fullfile(fileparts(mfilename('fullpath')),'_storage','ALL'))
pptfolder = fullfile(fileparts(mfilename('fullpath')),'report');

doexp0 = true;
doexp1 = true;
doexp2 = true;
doexp3 = true;
doexp5 = true;

if doexp0
    
    %[ppt0,op0]=openppt(fullfile(pptfolder,'plots_exp_0'),true);
    op0 = nan;
    
    A = exp0(op0,'canal',100,0,'pretimed');
    B = exp0(op0,'canal',100,0,'gaussian');
    C = exp0(op0,'canal',100,0.1,'pretimed');
    D = exp0(op0,'canal',100,0.1,'gaussian');
    
    % compute pretimed bandwidth for D
    D.setWindowType('pretimed');
    [total_band,band_o,band_i]=D.compute_total_bandwidth()
    
    % test half integer
    [to,ti]=A.segment_travel_times();
    
    abs_halfint_i = round([A.intersection.absoffseti]/50)*50;
    abs_halfint_o = round([A.intersection.absoffseto]/50)*50;
    
    [rel_halfint_o,rel_halfint_i] = A.absolute2reltaive(abs_halfint_o,abs_halfint_i,to,ti);
    
    for i=1:length(A.intersection)
       A.intersection(i).absoffseto = abs_halfint_o(i);
       A.intersection(i).absoffseti = abs_halfint_i(i);
       A.intersection(i).reloffseto = rel_halfint_o(i);
       A.intersection(i).reloffseti = rel_halfint_i(i);
    end
    
    figure('Position',[403 50 560 616],'visible','on')
    A.plot();
    addslide(op0,'half integer optimization')
    close
    
    %[total_band,band_o,band_i] = A.compute_total_bandwidth();
     
    closeppt(ppt0,op0)
end

if doexp1
    
    [ppt1p,op1p]=openppt(fullfile(pptfolder,'plots_exp_1_pretimed'),true);
    [ppt1g,op1g]=openppt(fullfile(pptfolder,'plots_exp_1_gaussian'),true);
    
    runnumber = 0;
    for average_segment_length = ALL.average_segment_length
        for cycle = ALL.cycle
            for green_split = ALL.green_split
                for max_abs_delta = ALL.max_abs_delta
                    for speed = ALL.speed
                        for numsamples = ALL.numsamples
                            runnumber = runnumber + 1;
                            disp(num2str(runnumber))
                            exp1name = fullfile(fileparts(mfilename('fullpath')),'_storage',...
                                ['exp1_'  ...
                                num2str(average_segment_length) '_' ...
                                num2str(cycle) '_' ...
                                num2str(green_split*100,'%.2d') '_' ...
                                num2str(max_abs_delta*100,'%.2d') '_' ...
                                num2str(speed) '_' ...
                                num2str(numsamples) ] );
                            
                            if(~exist([exp1name '.mat'],'file'))
                                error(['didn''t find ' exp1name])
                            end
                            
                            load(exp1name)
                            
                            tit = ['seg_length=' num2str(average_segment_length) ' ' ...
                                'cycle=' num2str(cycle) ' ' ...
                                'green_split=' num2str(green_split*100,'%.2d') ' ' ...
                                'delta_max=' num2str(max_abs_delta*100,'%.2d') ' ' ...
                                'speed=' num2str(speed) ' ' ...
                                'numsamples=' num2str(numsamples) ];
                            
                            make_exp1_figure(pretimed_ip.totalband,num_intersections,numsamples)
                            addslide(op1p,tit)
                            close
                            
                            make_exp1_figure(gaussian.totalband,num_intersections,numsamples)
                            addslide(op1g,tit)
                            close
                            
                        end
                    end
                end
            end
        end
    end
    closeppt(ppt1p,op1p)
    closeppt(ppt1g,op1g)
    
end

if doexp2
        
    runnumber = 0;
    for average_segment_length = ALL.average_segment_length
        for cycle = ALL.cycle
            for green_split = ALL.green_split
                for max_abs_delta = ALL.max_abs_delta
                    for speed = ALL.speed
                        for numsamples = ALL.numsamples
                            runnumber = runnumber + 1;
                            disp(num2str(runnumber))
                            
                            exp1name = fullfile(fileparts(mfilename('fullpath')),'_storage',...
                                ['exp1_'  ...
                                num2str(average_segment_length) '_' ...
                                num2str(cycle) '_' ...
                                num2str(green_split*100,'%.2d') '_' ...
                                num2str(max_abs_delta*100,'%.2d') '_' ...
                                num2str(speed) '_' ...
                                num2str(numsamples) ] );
                            
                            exp2(exp1name)
                        end
                    end
                end
            end
        end
    end
end

if doexp3
    
    [ppt3p,op3p]=openppt(fullfile(pptfolder,'plots_exp_3_pretimed'),true);
    % [ppt3g,op3g]=openppt(fullfile('plots_exp_3_gaussian'),true);
    % [ppt3z,op3z]=openppt(fullfile('plots_exp_3_compare'),true);
    
    runnumber = 0;
    for average_segment_length = ALL.average_segment_length
        for cycle = ALL.cycle
            for green_split = ALL.green_split
                for max_abs_delta = ALL.max_abs_delta
                    for speed = ALL.speed
                        for numsamples = ALL.numsamples
                            runnumber = runnumber+1;
                            disp(['********** RUN NUMBER ' num2str(runnumber) ' ******************'])
                            
                            exp1name = fullfile(fileparts(mfilename('fullpath')),'_storage',...
                                ['exp1_'  ...
                                num2str(average_segment_length) '_' ...
                                num2str(cycle) '_' ...
                                num2str(green_split*100,'%.2d') '_' ...
                                num2str(max_abs_delta*100,'%.2d') '_' ...
                                num2str(speed) '_' ...
                                num2str(numsamples)] );
                            
                            exp3name = regexprep(exp1name,'exp1','exp3');
                            
                            if(~exist([exp1name '.mat'],'file') | ~exist([exp3name '.mat'],'file'))
                                error('run exp1 and exp3 first')
                            end
                            
                            load(exp1name)
                            load(exp3name)
                            
                            tit = ['seg_length=' num2str(average_segment_length) ' ' ...
                                'cycle=' num2str(cycle) ' ' ...
                                'green_split=' num2str(green_split*100,'%.2d') ' ' ...
                                'delta_max=' num2str(max_abs_delta*100,'%.2d') ' ' ...
                                'speed=' num2str(speed) ' ' ...
                                'numsamples=' num2str(numsamples)];
                            
                            figure('Position',[180   340   560   326],'Visible','off')
                            n = length(num_intersections);
                            plot(1:n,mean(pert_pretimed_simplex,2),'k-','LineWidth',2), hold on
                            plot(1:n,mean(pert_pretimed_ip,2),'k:','LineWidth',2)
                            %     values=reshape(pretimed_simplex_var',1,n*numsamples);
                            %     tags = reshape( ones(numsamples,1)*num_intersections , 1, n*numsamples);
                            %     h=boxplot(values,tags,'notch','on','symbol','ko');
                            %     set(h(1:6,:),'LineWidth',2.5,'Color','k')
                            grid
                            set(gca,'xtick',1:n)
                            set(gca,'xticklabel',{num2str(num_intersections')})
                            set(gca,'FontSize',12)
%                             xlabel('number of intersections')
%                             ylabel('\Delta bandwidth [seconds]')
                            legend('simplex','ip')
                            
                            addslide(op3p,tit)
                            close
                            
                            %                         figure('Position',[180   246   896   420],'Visible','on')
                            %                         n = length(num_intersections);
                            %                         values=reshape(pert_gaussian',1,n*numsamples);
                            %                         tags = reshape( ones(numsamples,1)*num_intersections , 1, n*numsamples);
                            %                         h=boxplot(values,tags,'notch','on','symbol','ko');
                            %                         set(h(1:6,:),'LineWidth',2.5,'Color','k')
                            %                         set(gca,'xtick',1:n)
                            %                         set(gca,'xticklabel',{num2str(num_intersections')})
                            %                         set(gca,'FontSize',12)
                            %                         grid
                            %                         xlabel('number of intersections')
                            %                         ylabel('\Delta bandwidth [seconds]')
                            %                         addslide(op3g,tit)
                            %                         close
                            
                            %                         figure('Position',[180 243 813 423],'Visible','off')
                            %                         plot(1:n,mean(pretimed_simplex_var,2),'Color',0*[1 1 1],'LineWidth',2,'LineStyle','-')
                            %                         hold on
                            %                         plot(1:n,mean(pretimed_ip_var,2),'Color','r','LineWidth',2,'LineStyle','-')
                            %                         plot(1:n,mean(gaussian_ip_var,2),'Color',0*[1 1 1],'LineWidth',2,'LineStyle','--')
                            %                         plot(1:n,mean(gaussian_simplex_var,2),'Color',0*[1 1 1],'LineWidth',2,'LineStyle','--')
                            %                         legend('simplex','ip','gaussian')
                            %                         grid
                            %                         set(gca,'xtick',1:n)
                            %                         set(gca,'xticklabel',{num2str(num_intersections')})
                            %                         set(gca,'FontSize',12)
                            %                         grid
                            %                         xlabel('number of intersections')
                            %                         ylabel('\Delta bandwidth [seconds]')
                            %                         grid
                            %                         addslide(op3z,tit)
                            %                         close
                        end
                    end
                end
            end
        end
    end
    
    closeppt(ppt3p,op3p)
    % closeppt(ppt3g,op3g)
    % closeppt(ppt3z,op3z)
    
end

if doexp5
    
    load(fullfile(fileparts(mfilename('fullpath')),'_storage','ALL'))
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
                    
                    exp1name = fullfile(fileparts(mfilename('fullpath')),'_storage',...
                        ['exp1_'  ...
                        num2str(average_segment_length) '_' ...
                        num2str(cycle) '_' ...
                        num2str(green_split*100,'%.2d') '_' ...
                        num2str(max_abs_delta*100,'%.2d') '_' ...
                        num2str(speed) '_' ...
                        num2str(numsamples) ]);

                    load(exp1name)
                    
                    PretimedBand(i1,i2,i3,i4,:) = mean(pretimed_ip.totalband,2);
                    GaussianBand(i1,i2,i3,i4,:) = mean(gaussian.totalband,2);
                end
            end
        end
    end
    clear pretimed_* gaussian
    
    for i=1:length(ALL.cycle)
        leg{i} = [num2str(ALL.cycle(i)) ' sec. cycle'];
    end
    
    [ppt5,op5]=openppt(fullfile(pptfolder,'plots_exp_5'),true);
    for i1 = 1:n1
        for i3 = 1:n3
            for i4 = 1:n4
                average_segment_length = ALL.average_segment_length(i1);
                green_split = ALL.green_split(i3);
                max_abs_delta = ALL.max_abs_delta(i4);
                
                figure
                h=plot(num_intersections,fliplr(reshape(PretimedBand(i1,:,i3,i4,:),n2,7)'),'LineWidth',2);
                legend(fliplr(leg))
                grid
                set(h(4),'Color',0.0*[1 1 1],'Marker','+','MarkerSize',10)
                set(h(3),'Color',0.2*[1 1 1],'Marker','o','MarkerSize',10)
                set(h(2),'Color',0.4*[1 1 1],'Marker','x','MarkerSize',12)
                set(h(1),'Color',0.6*[1 1 1],'Marker','d','MarkerSize',10)
                xlabel('number of intersections')
                ylabel('total bandwidth [seconds]')
                ylim = get(gca,'YLim');
                set(gca,'YLim',[0.9*ylim(1) ylim(2)])
                set(gca,'XLim',[1.5 8.5])
                
                addslide(op5,['length=' num2str(average_segment_length) ' green_split=' num2str(green_split) ' delta=' num2str(max_abs_delta)])
                close
            end
            
        end
    end
    closeppt(ppt5,op5) 
end

disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=make_exp1_figure(totalband,num_intersections,numsamples)

figure('Position',[180   246   896   420],'Visible','off')
n = length(num_intersections);
plot(1:n,median(totalband,2),'Color',0.5*[1 1 1],'LineWidth',3,'LineStyle',':')
hold on
values=reshape(totalband',1,n*numsamples);
tags = reshape( ones(numsamples,1)*num_intersections , 1, n*numsamples);
h=boxplot(values,tags,'notch','on','symbol','ko');
set(h(1:6,:),'LineWidth',2.5,'Color','k')
set(gca,'xtick',1:n)
set(gca,'xticklabel',{num2str(num_intersections')})
set(gca,'FontSize',14)
grid
% xlabel('number of intersections')
% ylabel('bandwidth [seconds]')
