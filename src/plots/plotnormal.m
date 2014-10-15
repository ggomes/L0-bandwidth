function []=plotnormal()
% isolated function for playing with normal distributions
% and testing the method for maximization of Gaussian mixture.

close all

n = 3;
average_green = 10;
C = 30;
numpoints = 100;
dx = 2*C/(numpoints-1);

% individual mean and variance
sigma_in  = average_green*rand(1,n);
sigma_out = average_green*rand(1,n);
mu_in  = [0,C*rand(1,n-1)-C/2];
mu_out = [0,C*rand(1,n-1)-C/2];
delta = C*rand(1,n)-C/2;
x = linspace(-C,C,numpoints);


% individual gaussians
N_out = zeros(n,numpoints);
N_in = zeros(n,numpoints);
for i=1:n
    N_out(i,:) = gaussian(mu_out(i),sigma_out(i),x);
    N_in(i,:) = gaussian(mu_in(i),sigma_in(i),x);
end

% evaluate product of gaussians
[p_in_gamma,p_in_mu,p_in_sigma]=evaluateproduct(mu_in,sigma_in);
[p_out_gamma,p_out_mu,p_out_sigma]=evaluateproduct(mu_out,sigma_out);


% plot and check formula
figure('Position',[139 39 560 629])
for i=1:n
    subplot(n+1,1,i)
    subplotf(n+1,i,'','',x,N_in(i,:),'','k','',2);
    vline(mu_in(i),'r')
    grid
end
% subplotf(n+1,n+1,'','',x,prod(N_in,1),'','k','',2);
subplotf(n+1,n+1,'','',x,p_in_gamma*gaussian(p_in_mu,p_in_sigma,x),'','r--','',2);
vline(p_in_mu)
grid

% plot gamma function with n=3
mu2 = x;
mu3 = x;
for i=1:length(mu2)
    for j=1:length(mu3)
        m_in = [0 mu2(i) mu3(j)];
        m_out = m_in+delta;
        m_out(1)=0;
        gamma_in(i,j)=evaluategamma(m_in,sigma_in);
        gamma_out(i,j)=evaluategamma(m_out,sigma_out);
    end
end

gamma = gamma_in+gamma_out;

% numerical gradient
% [U_in,V_in] = gradient(gamma_in,dx);
% [U_out,V_out] = gradient(gamma_out,dx);
% [U,V] = gradient(gamma,dx);

% analytical gradient
for i=1:length(mu2)
    for j=1:length(mu3)
        m_in = [0 mu2(i) mu3(j)];
        d = evaluateDgamma(m_in,sigma_in);
        Va_in(i,j) = d(1);
        Ua_in(i,j) = d(2);
        
        m_out = m_in+delta;
        m_out(1)=0;
        d = evaluateDgamma(m_out,sigma_out);
        Va_out(i,j) = d(1);
        Ua_out(i,j) = d(2);
    end
end

Ua = Ua_in+Ua_out;
Va = Va_in+Va_out;

% plot gamma function
figure('Position',[95 260  1265 420])

subplot(121)
contour(mu2,mu3,gamma_in,'LineWidth',2,'Color',[0 0 0]);
hold on
contour(mu2,mu3,gamma_out,'LineWidth',2,'Color',[0 0 0])
quiver(mu2,mu3,Ua_in,Va_in)
quiver(mu2,mu3,Ua_out,Va_out)

% sum gamma
subplot(122)
contour(mu2,mu3,gamma,'LineWidth',2,'Color',[0 0 0]);
hold on
quiver(mu2,mu3,Ua,Va)

% initial condition is highest peak
gammapeak_in = evaluategamma([0 0 0],sigma_in);
gammapeak_out = evaluategamma([0 0 0],sigma_out);
if(gammapeak_in>gammapeak_out)
    p = [0 0];
else
    p = -delta(2:3);
end

% epsilon = 0.1;
% for i=1:10
%     g = 
% end

subplot(121), hold on
plot(ic(2),ic(1),'.','MarkerSize',20,'Color',[1 0 0])
subplot(122), hold on
plot(ic(2),ic(1),'.','MarkerSize',20,'Color',[1 0 0])

% -------------------------------------------------------------------------
function [y]=gaussian(mu,sigma,x)
y=exp(-(x-mu).^2/2/sigma/sigma)/sqrt(2*pi)/sigma;

% -------------------------------------------------------------------------
function [p_gamma,p_mu,p_sigma]=evaluateproduct(mu,sigma)
% evaluate product formula

sumsigmainvsquare = sum(1./(sigma.^2));
p_gamma = evaluategamma(mu,sigma);
p_mu = sum( mu./(sigma.^2) )/sumsigmainvsquare;
p_sigma = sqrt(1/sumsigmainvsquare);

% -------------------------------------------------------------------------
function [gamma]=evaluategamma(mu,sigma)

n = length(mu);
s = 0;
prodsigmasquare = prod(sigma.^2);
sumsigmainvsquare = sum(1./(sigma.^2));
for i=1:n-1
   for j=i+1:n 
        s = s+((mu(i)-mu(j))/sigma(i)/sigma(j))^2;
   end
end
gamma = exp(-s/sumsigmainvsquare/2)/(2*pi)^((n-1)/2)/sqrt(prodsigmasquare*sumsigmainvsquare);

% -------------------------------------------------------------------------
function [dgamma]=evaluateDgamma(mu,sigma)

sigmabarsquare = sum(1./(sigma.^2));
mubar = sum(mu./(sigma.^2))/sigmabarsquare;
gamma = evaluategamma(mu,sigma);
dgamma = gamma*(mubar-mu)./(sigma.^2);
dgamma = dgamma(2:end);

