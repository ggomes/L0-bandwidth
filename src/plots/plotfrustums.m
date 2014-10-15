function []=plotfrustums()
% isolated function for playing with 3D frustums

W = 20;
N = 150;

w2 = ones(N,1)*linspace(-W,W,N);
w3 = linspace(-W,W,N)'*ones(1,N);
A= evalfrustum(5,6,8,0,0,1);

for alpha = 1
    
    B= evalfrustum(5,6,8,0,-8,alpha);
    clf
    h=surf(A+B);
    set(h,'EdgeAlpha',0)
    view(90,-90)
    pause(0.1)
end

    function [Z]=evalfrustum(g1,g2,g3,dw2,dw3,alpha)
        g12 = alpha*(g1+g2)/2;
        g13 = alpha*(g1+g3)/2;
        g23 = alpha*(g2+g3)/2;
        mw2 = alpha*(w2+dw2);
        mw3 = alpha*(w3+dw3);
        X={mw2+g12,-mw2+g12,mw3+g13,-mw3+g13,mw2-mw3+g23,mw3-mw2+g23,g1,g2,g3};
        Z=X{1};
        for ii=2:length(X)
            Z=min(Z,X{ii});
        end
        Z=max(Z,0);
    end

end
