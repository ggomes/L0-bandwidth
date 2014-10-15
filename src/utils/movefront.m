function []=movefront(h)

a = get(gca,'Children');
p = find(a==h);
v = a;
v(p) = [];
v = [h;v];
set(gca,'Children',v)