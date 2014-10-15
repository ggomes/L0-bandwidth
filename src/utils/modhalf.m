function [b]=modhalf(a,cycle)

b = mod( a + cycle/2 , cycle ) - cycle/2;
