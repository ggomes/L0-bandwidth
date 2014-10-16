function [y]=eval_gauss(x,mu,sigma,gamma)
y = gamma * exp( -((x-mu).^2)/2/sigma^2 )/sqrt(2*pi)/sigma;
