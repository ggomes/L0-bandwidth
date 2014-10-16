function [Mu,Sigma] = gaussian_product(mu,sigma)

SigmaSquare = 1/sum(sigma.^-2);               % Eq (73)
Mu = SigmaSquare * sum( mu .* sigma.^-2 );    % Eq (72)
Sigma = sqrt(SigmaSquare);
