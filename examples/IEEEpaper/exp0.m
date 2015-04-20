function [A]=exp0(op,facility,cycle,delta,windowtype)
% Experiment #0
% -------------
% solve for single facility

if(strcmp(facility,'canal'))
    A = load_Canal(cycle,delta,windowtype);
end
if(strcmp(facility,'sanpablo'))
    A = load_SanPablo(cycle,delta,windowtype);
    A.remove_intersection('Allston');
    %A.remove_intersection('Grayson');
    A.remove_intersection('Dwight');
    A.remove_intersection('Addison');
end

figure('Position',[403 50 560 616],'visible','on')
A.optimize();
A.plot();

if(~isnan(op))
    addslide(op,[facility ' cycle=' num2str(cycle) ' delta=' num2str(delta)])
end


