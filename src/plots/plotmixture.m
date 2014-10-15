function [f]=plotmixture(bo,alphao,bi,alphai,eo,ei)

x = linspace(0,1);
yo = bo*exp(-alphao*x.^2);
yi = bi*exp(-alphai*(1-x).^2);
Y = yi+yo;
lambdao = 1/sqrt(2*alphao);
lambdai = 1/sqrt(2*alphai);
maxind = find(diff(double(diff(Y)>0))==-1)+1;

Y_eo = bo*exp(-alphao*eo.^2) + bi*exp(-alphai*(1-eo).^2);
Y_ei = bo*exp(-alphao*ei.^2) + bi*exp(-alphai*(1-ei).^2);

f=figure('Position',[223 336 740 330]);
jbfill(gcf,x,0*yo,yo,0.5*ones(1,3),0.5*ones(1,3),true,0.3);
jbfill(gcf,x,0*yi,yi,0.5*ones(1,3),0.5*ones(1,3),true,0.3);
hold on
set(gca,'YLim',[0,1.2*max(Y)])
h(1) = vline(lambdao);
h(2) = vline(1-lambdai);
set(h,'LineWidth',1,'Color',zeros(1,3),'LineStyle','--')
plot(x,Y,'k.-','LineWidth',1.5)
plot(x(maxind),Y(maxind),'.','MarkerSize',25,'LineWidth',2,'Color',0.5*ones(1,3))

plot(eo,Y_eo,'r.','MarkerSize',15,'LineWidth',2)
plot(eo(end),Y_eo(end),'mo','MarkerSize',8,'LineWidth',2)

plot(ei,Y_ei,'r.','MarkerSize',15,'LineWidth',2)
plot(ei(end),Y_ei(end),'mo','MarkerSize',8,'LineWidth',2)

set(gca,'XLim',[0 1])
