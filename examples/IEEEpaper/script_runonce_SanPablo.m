clear 
close all     

cycle = 100;
windowtype = 'gaussian';

A = load_SanPablo(cycle,windowtype);

A.remove_intersection('Allston');
%A.remove_intersection('Grayson');
A.remove_intersection('Dwight');
A.remove_intersection('Addison');

clf

A.optimize();
A.plot(gcf);

disp('done')




