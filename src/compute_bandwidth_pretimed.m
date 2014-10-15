function [band]=compute_bandwidth_pretimed(reloffset,green)

interval = reloffset*[1 1] + green/2*[-1 1];
a = max(interval(:,1));
b = min(interval(:,2));
if(b>a)
    band = b-a;
else
    band = 0;
end
