function plotincifield(hinci,event,nx,ny,nz,x,y,z,matxyincifield,matxyincifieldx,matxyincifieldy,matxyincifieldz,nprint)
  val = get(hinci,'Value');


switch val

case 2

figure(10)
set(10,'DefaultAxesFontName','Times')
set(10,'DefaultAxesFontSize',12)
set(10,'DefaultAxesFontWeight','Bold')
set(10,'DefaultTextfontName','Times')
set(10,'DefaultTextfontSize',12)
set(10,'DefaultTextfontWeight','Bold')
set(10,'Position',[0 0 1000 600])

uicontrol('Style', 'text', 'String', 'Choice of z','Position', [250 10 90 20]);
uicontrol('Style', 'text', 'String', 'm','Position', [640 10 15 20]);
uicontrol('Style', 'text', 'String', num2str(z(1),'%+1.2e\n'),'Position', [560 10 80 20]);

  subplot('position',[0.075 0.63 0.33 0.33])
  imagesc(x,y,matxyincifield(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('modulus')  
colorbar('vert')

  
  subplot('position',[0.65 0.63 0.33 0.33])
  imagesc(x,y,abs(matxyincifieldx(:,:,1)'));
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title(' x')  
colorbar('vert')
  
 subplot('position',[0.075 0.15 0.33 0.33])
  imagesc(x,y,abs(matxyincifieldy(:,:,1)'));
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title(' y')  
colorbar('vert')

subplot('position',[0.65 0.15 0.33 0.33])
imagesc(x,y,abs(matxyincifieldz(:,:,1)'));
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title(' z')  
colorbar('vert')


uicontrol('Style', 'slider', 'Min',1,'Max', nz,...
	  'val',1,'sliderstep',[1/(nz-1) 2/(nz-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotincifieldz,x,y,z,matxyincifield,matxyincifieldx,matxyincifieldy,matxyincifieldz,nprint});

case 3


figure(10)
set(10,'DefaultAxesFontName','Times')
set(10,'DefaultAxesFontSize',12)
set(10,'DefaultAxesFontWeight','Bold')
set(10,'DefaultTextfontName','Times')
set(10,'DefaultTextfontSize',12)
set(10,'DefaultTextfontWeight','Bold')
set(10,'Position',[0 0 1000 600])

uicontrol('Style', 'text', 'String', 'Choice of y','Position', [250 10 90 20]);
uicontrol('Style', 'text', 'String', 'm','Position', [640 10 15 20]);
uicontrol('Style', 'text', 'String', num2str(y(1),'%+1.2e\n'),'Position', [560 10 80 20]);


  subplot('position',[0.075 0.63 0.33 0.33])
  imagesc(x,z,permute(matxyincifield(:,1,:),[1 3 2])');
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title('modulus')  
colorbar('vert')

  
  subplot('position',[0.65 0.63 0.33 0.33])
  imagesc(x,z,abs(permute(matxyincifieldx(:,1,:),[1 3 2])'));
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title(' x')  
colorbar('vert')
  
 subplot('position',[0.075 0.15 0.33 0.33])
  imagesc(x,z,abs(permute(matxyincifieldy(:,1,:),[1 3 2])'));
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title(' y')  
colorbar('vert')


  subplot('position',[0.65 0.15 0.33 0.33])
  imagesc(x,z,abs(permute(matxyincifieldz(:,1,:),[1 3 2])'));
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title(' z')  
colorbar('vert')

uicontrol('Style', 'slider', 'Min',1,'Max', ny,...
	  'val',1,'sliderstep',[1/(ny-1) 2/(ny-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotincifieldy,x,y,z,matxyincifield,matxyincifieldx,matxyincifieldy,matxyincifieldz,nprint});


case 4


figure(10)
set(10,'DefaultAxesFontName','Times')
set(10,'DefaultAxesFontSize',12)
set(10,'DefaultAxesFontWeight','Bold')
set(10,'DefaultTextfontName','Times')
set(10,'DefaultTextfontSize',12)
set(10,'DefaultTextfontWeight','Bold')
set(10,'Position',[0 0 1000 600])

uicontrol('Style', 'text', 'String', 'Choice of x','Position', [250 10 90 20]);
uicontrol('Style', 'text', 'String', 'm','Position', [640 10 15 20]);
uicontrol('Style', 'text', 'String', num2str(x(1),'%+1.2e\n'),'Position', [560 10 80 20]);

  subplot('position',[0.075 0.63 0.33 0.33])
  imagesc(y,z,permute(matxyincifield(1,:,:),[2 3 1])');
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title('modulus')  
colorbar('vert')

  
  subplot('position',[0.65 0.63 0.33 0.33])
  imagesc(y,z,abs(permute(matxyincifieldx(1,:,:),[2 3 1])'));
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title(' x')  
colorbar('vert')
  
 subplot('position',[0.075 0.15 0.33 0.33])
  imagesc(y,z,abs(permute(matxyincifieldy(1,:,:),[2 3 1])'));
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title('y')  
colorbar('vert')


  subplot('position',[0.65 0.15 0.33 0.33])
  imagesc(y,z,abs(permute(matxyincifieldz(1,:,:),[2 3 1])'));
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title(' z')  
colorbar('vert')

  

uicontrol('Style', 'slider', 'Min',1,'Max', nx,...
	  'val',1,'sliderstep',[1/(nx-1) 2/(nx-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotincifieldx,x,y,z,matxyincifield,matxyincifieldx,matxyincifieldy,matxyincifieldz,nprint});


end;

if (nprint == 1)
print('-f10','incident','-depsc')
end
