function [ ] = paper_plots(  )

add_dependencies()

plot_scalar_mixture

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function []=plot_scalar_mixture()

gamma0= 1;
gamma1= 1;

lambda0 = 0.3;
lambda1 = 0.5;

epsilon = linspace(0,1);

xi0 = (gamma0/sqrt(2*pi)/lambda0)*exp(-(1/2/lambda0^2)*epsilon.^2);
xi1 = (gamma1/sqrt(2*pi)/lambda1)*exp(-(1/2/lambda1^2)*(1-epsilon).^2);
xi = xi0+xi1;

% find max
maxind = find(diff(double(diff(xi)>0))==-1)+1;

figure('Position',[223   336   740   330])
jbfill(gcf,epsilon,0*xi0,xi0,0.5*ones(1,3),0.5*ones(1,3),true,0.3);
jbfill(gcf,epsilon,0*xi1,xi1,0.5*ones(1,3),0.5*ones(1,3),true,0.3);
hold on
set(gca,'YLim',[0,1.2*max(xi)])
h(1) = vline(lambda0);
h(2) = vline(1-lambda1);
set(h,'LineWidth',1,'Color',zeros(1,3),'LineStyle','--')
plot(epsilon,xi,'k','LineWidth',1.5)
plot(epsilon(maxind),xi(maxind),'.','MarkerSize',25,'LineWidth',2,'Color',0.5*ones(1,3))
set(gca,'XLim',[0 1])
