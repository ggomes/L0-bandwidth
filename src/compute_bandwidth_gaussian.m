function [band]=compute_bandwidth_gaussian(reloffset,sigma,gamma)

n = length(sigma);
A = zeros(n-1);
for i=2:n
    for j=2:n
        A(i-1,j-1) = (sigma(i)*sigma(j))^-2;
    end
end
band_sigma = diag(sigma(2:end).^-2) - A/sum(sigma.^-2);
band_b = prod(gamma)/sqrt( ((2*pi)^(n-1)) * sum(sigma.^-2) * prod(sigma.^2) );

reloffset_sub = reloffset(2:end);
band = band_b*exp(-reloffset_sub'*band_sigma*reloffset_sub/2);
