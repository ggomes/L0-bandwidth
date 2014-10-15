function [Sigma,bo] = integral_gaussian_product(gamma,sigma)

n = length(gamma);
A = zeros(n-1);
for i=2:n
    for j=2:n
        A(i-1,j-1) = (sigma(i)*sigma(j))^-2;
    end
end
Sigma = diag(sigma(2:end).^-2) - A/sum(sigma.^-2);

bo = prod(gamma)/sqrt( ((2*pi)^(n-1)) * sum(sigma.^-2) * prod(sigma.^2) );
