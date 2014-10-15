function [h3]=subplotf(n,p,ylab,xlab,xdata,ydata,tit,c,leg,x)
% (n,p,ylab,xlab,xdata,ydata,tit,c,leg,x)

labelfontsize = 12;
titlefontsize = 12;
legendfontsize = 10;

h= 2*(43-n)/n ;

q=n-p;  % row number counted from bottom
u=8+(h+2)*q;

% hold on
h1=subplot('position',[0.1 u/100 0.88 h/100]);

if(nargin==2)
    h3=h1;
    return
end

hold on
if(isempty(c))
    h3=plot(xdata,ydata);
else
    h3=plot(xdata,ydata,c);
end
h2=ylabel(ylab,'FontSize',labelfontsize);
set(h2,'VerticalAlignment','bottom')

if(q>0)
    set(h1,'XTickLabel',[])
else
    h=xlabel(xlab);
    set(h,'FontSize',labelfontsize)
    if(~isempty(leg))
        h=legend(leg,1);
        set(h,'FontSize',legendfontsize)
    end
end

if((p==1))
    if(~isempty(tit))
        title(tit,'FontSize',titlefontsize,'FontWeight','bold')
    end
end

set(h3,'LineWidth',x)