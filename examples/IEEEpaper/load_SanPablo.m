function [A]=load_SanPablo(cycle,delta,windowtype)

A = class_artery(cycle,windowtype);

% I = class_intersection('Ashby',windowtype,35/110,31/110,-21/110);    % (gbar,g,d)
I = class_intersection('Ashby',windowtype,35/110,31/110,delta);            % (gbar,g,d)
A.add_intersection(I);
A.add_segment(1580,30,30);

I = class_intersection('Grayson',windowtype,71/110,71/110,delta);
A.add_intersection(I);
A.add_segment(1875,30,30);

I = class_intersection('Dwight',windowtype,37/110,37/110,delta);
A.add_intersection(I);
A.add_segment(2035,30,30);

I = class_intersection('Allston',windowtype,25/60,25/60,delta);
A.add_intersection(I);
A.add_segment(585,30,30);

I = class_intersection('Addison',windowtype,26/60,26/60,delta);
A.add_intersection(I);
A.add_segment(500,30,30);

I = class_intersection('University',windowtype,43/120,43/120,delta);
A.add_intersection(I);
            
A.initialize();