clf
n = 100;
sigma = 0.05;
gamma = 1;
verts = [ [ zeros(n,1) linspace(0,1,n)' ] ; ...
          [  ones(n,1) linspace(0,1,n)' ] ];
verts(:,[2 1]) = verts(:,[1 2]);
ind = verts(:,2)==1;
verts(ind,1)=verts(ind,1)+3;
faces = [[1:n-1]' [2:n]' [n+2:2*n]' [n+1:2*n-1]'];
     
p = patch('Faces',faces,'Vertices',verts,'EdgeColor','none');

clear cdata
set(gca,'CLim',[0 1])
cdata = linspace(0,1,n-1)';
cdata = gamma*(1-exp( -(cdata-0.5).^2/sigma ));
set(p,'FaceColor','flat','FaceVertexCData',cdata,'CDataMapping','scaled')

colormap('gray')