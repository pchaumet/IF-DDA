function plotforce(hlocal,event,nx,ny,nz,x,y,z,xx,yy,zz,matxyforcex,matxyforcey,matxyforcez,nprint)

val = get(hlocal,'Value');

switch val

case 2

figure(200)
set(200,'DefaultAxesFontName','Times')
set(200,'DefaultAxesFontSize',12)
set(200,'DefaultAxesFontWeight','Bold')
set(200,'DefaultTextfontName','Times')
set(200,'DefaultTextfontSize',12)
set(200,'DefaultTextfontWeight','Bold')
set(200,'Position',[0 0 1000 600])

uicontrol('Style', 'text', 'String', 'Choice of z','Position', [250 15 90 20]);
uicontrol('Style', 'text', 'String', 'm','Position', [640 15 15 20]);
uicontrol('Style', 'text', 'String', num2str(z(1),'%+1.2e\n'),'Position', [560 15 80 20]);

scale=1
subplot('position',[0.2 0.15 0.6 0.7])

quiver(xx(:,:,1),yy(:,:,1),matxyforcex(:,:,1),matxyforcey(:,:,1),scale)

axis equal  
xlabel('x')
ylabel('y')

uicontrol('Style', 'slider', 'Min',1,'Max', nz,...
'val',1,'sliderstep',[1/(nz-1) 2/(nz-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotforcez,xx,yy,z,matxyforcex,matxyforcey,nprint});


case 3


figure(200)
set(200,'DefaultAxesFontName','Times')
set(200,'DefaultAxesFontSize',12)
set(200,'DefaultAxesFontWeight','Bold')
set(200,'DefaultTextfontName','Times')
set(200,'DefaultTextfontSize',12)
set(200,'DefaultTextfontWeight','Bold')
set(200,'Position',[0 0 1000 600])

uicontrol('Style', 'text', 'String', 'Choice of y','Position', [250 15 90 20]);
uicontrol('Style', 'text', 'String', 'm','Position', [640 15 15 20]);
uicontrol('Style', 'text', 'String', num2str(y(1),'%+1.2e\n'),'Position', [560 15 80 20]);

scale=1
subplot('position',[0.2 0.15 0.6 0.7])

  quiver(xx(:,1,:),zz(:,1,:),matxyforcex(:,1,:),matxyforcez(:,1,:),scale)

axis equal  
xlabel('x')
ylabel('z')

uicontrol('Style', 'slider', 'Min',1,'Max', ny,...
'val',1,'sliderstep',[1/(ny-1) 2/(ny-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotforcey,xx,zz,y,matxyforcex,matxyforcez,nprint});



case 4


figure(200)
set(200,'DefaultAxesFontName','Times')
set(200,'DefaultAxesFontSize',12)
set(200,'DefaultAxesFontWeight','Bold')
set(200,'DefaultTextfontName','Times')
set(200,'DefaultTextfontSize',12)
set(200,'DefaultTextfontWeight','Bold')
set(200,'Position',[0 0 1000 600])

uicontrol('Style', 'text', 'String', 'Choice of x','Position', [250 15 90 20]);
uicontrol('Style', 'text', 'String', 'm','Position', [640 15 15 20]);
uicontrol('Style', 'text', 'String', num2str(x(1),'%+1.2e\n'),'Position', [560 15 80 20]);

scale=1
subplot('position',[0.2 0.15 0.6 0.7])
  
  quiver(yy(1,:,:),zz(1,:,:),matxyforcey(1,:,:),matxyforcez(1,:,:),scale)

axis equal  
xlabel('y')
ylabel('z')

uicontrol('Style', 'slider', 'Min',1,'Max', nx,...
'val',1,'sliderstep',[1/(nx-1) 2/(nx-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotforcex,yy,zz,x,matxyforcey,matxyforcez,nprint});



end;

if (nprint == 1)
print('-f200','force2d','-depsc')
end
