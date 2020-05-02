function plotmacrofield(hmacro,event,nx,ny,nz,x,y,z,matxymacrofield,matxymacrofieldx,matxymacrofieldy,matxymacrofieldz,nprint)
  val = get(hmacro,'Value');


switch val

case 2

figure(30)
set(30,'DefaultAxesFontName','Times')
set(30,'DefaultAxesFontSize',12)
set(30,'DefaultAxesFontWeight','Bold')
set(30,'DefaultTextfontName','Times')
set(30,'DefaultTextfontSize',12)
set(30,'DefaultTextfontWeight','Bold')
set(30,'Position',[0 0 1000 600])

uicontrol('Style', 'text', 'String', 'Choice of z','Position', [250 10 90 20]);
uicontrol('Style', 'text', 'String', 'm','Position', [640 10 15 20]);
uicontrol('Style', 'text', 'String', num2str(z(1),'%+1.2e\n'),'Position', [560 10 80 20]);

  subplot('position',[0.075 0.63 0.33 0.33])
  imagesc(x,y,matxymacrofield(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('modulus')  
colorbar('vert')

  
  subplot('position',[0.65 0.63 0.33 0.33])
  imagesc(x,y,abs(matxymacrofieldx(:,:,1)'));
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title(' x')  
colorbar('vert')
  
 subplot('position',[0.075 0.15 0.33 0.33])
  imagesc(x,y,abs(matxymacrofieldy(:,:,1)'));
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title(' y')  
colorbar('vert')

%text(-0.4,0.8,'Macro field','units','Normalized','rotation',90)


  subplot('position',[0.65 0.15 0.33 0.33])
  imagesc(x,y,abs(matxymacrofieldz(:,:,1)'));
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title(' z')  
colorbar('vert')


uicontrol('Style', 'slider', 'Min',1,'Max', nz,...
	  'val',1,'sliderstep',[1/(nz-1) 2/(nz-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotmacrofieldz,x,y,z,matxymacrofield,matxymacrofieldx,matxymacrofieldy,matxymacrofieldz,nprint});

case 3


figure(30)
set(30,'DefaultAxesFontName','Times')
set(30,'DefaultAxesFontSize',12)
set(30,'DefaultAxesFontWeight','Bold')
set(30,'DefaultTextfontName','Times')
set(30,'DefaultTextfontSize',12)
set(30,'DefaultTextfontWeight','Bold')
set(30,'Position',[0 0 1000 600])

uicontrol('Style', 'text', 'String', 'Choice of y','Position', [250 10 90 20]);
uicontrol('Style', 'text', 'String', 'm','Position', [640 10 15 20]);
uicontrol('Style', 'text', 'String', num2str(y(1),'%+1.2e\n'),'Position', [560 10 80 20]);


  subplot('position',[0.075 0.63 0.33 0.33])
  imagesc(x,z,permute(matxymacrofield(:,1,:),[1 3 2])');
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title('modulus')  
colorbar('vert')

  
  subplot('position',[0.65 0.63 0.33 0.33])
  imagesc(x,z,abs(permute(matxymacrofieldx(:,1,:),[1 3 2])'));
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title(' x')  
colorbar('vert')
  
 subplot('position',[0.075 0.15 0.33 0.33])
  imagesc(x,z,abs(permute(matxymacrofieldy(:,1,:),[1 3 2])'));
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title(' y')  
colorbar('vert')


  subplot('position',[0.65 0.15 0.33 0.33])
  imagesc(x,z,abs(permute(matxymacrofieldz(:,1,:),[1 3 2])'));
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title(' z')  
colorbar('vert')

uicontrol('Style', 'slider', 'Min',1,'Max', ny,...
	  'val',1,'sliderstep',[1/(ny-1) 2/(ny-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotmacrofieldy,x,y,z,matxymacrofield,matxymacrofieldx,matxymacrofieldy,matxymacrofieldz,nprint});


case 4


figure(30)
set(30,'DefaultAxesFontName','Times')
set(30,'DefaultAxesFontSize',12)
set(30,'DefaultAxesFontWeight','Bold')
set(30,'DefaultTextfontName','Times')
set(30,'DefaultTextfontSize',12)
set(30,'DefaultTextfontWeight','Bold')
set(30,'Position',[0 0 1000 600])

uicontrol('Style', 'text', 'String', 'Choice of x','Position', [250 10 90 20]);
uicontrol('Style', 'text', 'String', 'm','Position', [640 10 15 20]);
uicontrol('Style', 'text', 'String', num2str(x(1),'%+1.2e\n'),'Position', [560 10 80 20]);

  subplot('position',[0.075 0.63 0.33 0.33])
  imagesc(y,z,permute(matxymacrofield(1,:,:),[2 3 1])');
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title('modulus')  
colorbar('vert')

  
  subplot('position',[0.65 0.63 0.33 0.33])
  imagesc(y,z,abs(permute(matxymacrofieldx(1,:,:),[2 3 1])'));
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title(' x')  
colorbar('vert')
  
 subplot('position',[0.075 0.15 0.33 0.33])
  imagesc(y,z,abs(permute(matxymacrofieldy(1,:,:),[2 3 1])'));
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title('y')  
colorbar('vert')


  subplot('position',[0.65 0.15 0.33 0.33])
  imagesc(y,z,abs(permute(matxymacrofieldz(1,:,:),[2 3 1])'));
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title(' z')  
colorbar('vert')

  

uicontrol('Style', 'slider', 'Min',1,'Max', nx,...
	  'val',1,'sliderstep',[1/(nx-1) 2/(nx-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotmacrofieldx,x,y,z,matxymacrofield,matxymacrofieldx,matxymacrofieldy,matxymacrofieldz,nprint});


end;

if (nprint == 1)
print('-f30','macroscopic','-depsc')
end
