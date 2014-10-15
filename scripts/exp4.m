function []=exp4(facility,delta,windowtype)
% Experiment #4
% -------------
% For gaussian windows, do the resulting offsets depend on the window
% sigmas?
% Plot offsets versus a global sigma multiplier.
% As it turns out, the looks a lot like a step function with the 
% limit as mult->0 selecting the largest of the two peaks, and 
% the limit as mult->inft having an analytical form. 

add_dependencies()

if(strcmp(facility,'canal'))
    cycle = 200;
    A = load_Canal(cycle,delta,windowtype);
end
if(strcmp(facility,'sanpablo'))
    cycle = 120;
    A = load_SanPablo(cycle,windowtype);
    A.remove_intersection('Allston');
    %A.remove_intersection('Grayson');
    A.remove_intersection('Dwight');
    A.remove_intersection('Addison');
end

num_intersection = length(A.intersection);

mult = linspace(0.2,cycle/2);

sigma_o = rand(1,num_intersection);
sigma_i = rand(1,num_intersection);

gamma_o = sqrt(2*pi)*sigma_o;
gamma_i = sqrt(2*pi)*sigma_i;

for m=1:length(mult)
    
    for i=1:length(A.intersection)
        A.intersection(i).sigma_o = mult(m)*sigma_o(i);
        A.intersection(i).sigma_i = mult(m)*sigma_i(i);
        A.intersection(i).gamma_o = mult(m)*gamma_o(i);
        A.intersection(i).gamma_i = mult(m)*gamma_i(i);
    end
    
    A.optimize();
    
    AbsOffsetO(m,:) = [A.intersection.absoffseto];
    AbsOffsetI(m,:) = [A.intersection.absoffseti];
    
end

% round to the nearest 0.2 second
AbsOffsetI = round(AbsOffsetI*10)/10;
AbsOffsetO = round(AbsOffsetO*10)/10;

AbsOffsetI = modhalf(AbsOffsetI,cycle);
AbsOffsetO = modhalf(AbsOffsetO,cycle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Sigma_o,bo_o] = integral_gaussian_product(gamma_o,sigma_o);
[Sigma_i,bo_i] = integral_gaussian_product(gamma_i,sigma_i);
[to,ti]=A.segment_travel_times();
delta0 = A.translated_internal_offsets(to,ti);
delta_bold = delta0(1)-delta0(2:end);
alpha_o = delta_bold*Sigma_o*delta_bold' / 2;
alpha_i = delta_bold*Sigma_i*delta_bold' / 2;


e_inf = bo_i*alpha_i/(bo_i*alpha_i+bo_o*alpha_o);

omegaO_star = e_inf * delta_bold;
omegaO = [0;omegaO_star'];
omegaI = omegaO + delta0' - delta0(1);
omegaO = modhalf(omegaO,cycle);
omegaI = modhalf(omegaI,cycle);
[AbsOffsetO_inf,AbsOffsetI_inf] = A.relative2absolute(omegaO,omegaI,to,ti);
%%%%%%%%%%%%%%%%%%%
            
figure('Position',[403 117 620 549])
subplotf(2,1,'offset [sec]','multiplier',mult,AbsOffsetI,'','','',2);
set(gca,'YLim',1.05*cycle/2*[-1 1])
hline(AbsOffsetI_inf)
subplotf(2,2,'offset [sec]','multiplier',mult,AbsOffsetO,'','','',2);
set(gca,'YLim',1.05*cycle/2*[-1 1])
hline(AbsOffsetO_inf)

disp('done')

