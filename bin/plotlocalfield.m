function plotlocalfield(hlocal,event,nx,ny,nz,x,y,z,matxylocalfield,matxylocalfieldx,matxylocalfieldy,matxylocalfieldz)
  val = get(hlocal,'Value');


switch val

case 2

figure(20)
set(20,'DefaultAxesFontName','Times')
set(20,'DefaultAxesFontSize',12)
set(20,'DefaultAxesFontWeight','Bold')
set(20,'DefaultTextfontName','Times')
set(20,'DefaultTextfontSize',12)
set(20,'DefaultTextfontWeight','Bold')
set(20,'Position',[0 0 1000 600])

uicontrol('Style', 'text', 'String', 'Choice of z','Position', [250 10 90 20]);
uicontrol('Style', 'text', 'String', 'm','Position', [640 10 15 20]);
uicontrol('Style', 'text', 'String', num2str(z(1),'%+1.2e\n'),'Position', [560 10 80 20]);


%maxlocalz=max(max(max(matxylocalfield)));
%minlocalz=min(min(min(matxylocalfield)));

  subplot('position',[0.075 0.63 0.33 0.33])
  imagesc(x,y,matxylocalfield(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('modulus')  
colorbar('vert')
%caxis([minlocalz maxlocalz])
  
  subplot('position',[0.65 0.63 0.33 0.33])
  imagesc(x,y,abs(matxylocalfieldx(:,:,1)'));
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title(' x')  
colorbar('vert')
  
 subplot('position',[0.075 0.15 0.33 0.33])
  imagesc(x,y,abs(matxylocalfieldy(:,:,1)'));
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title(' y')  
colorbar('vert')

%text(-0.4,0.8,'Local field','units','Normalized','rotation',90)


  subplot('position',[0.65 0.15 0.33 0.33])
  imagesc(x,y,abs(matxylocalfieldz(:,:,1)'));
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title(' z')  
colorbar('vert')


uicontrol('Style', 'slider', 'Min',1,'Max', nz,...
	  'val',1,'sliderstep',[1/(nz-1) 2/(nz-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotlocalfieldz,x,y,z,matxylocalfield,matxylocalfieldx,matxylocalfieldy,matxylocalfieldz});

case 3


figure(20)
set(20,'DefaultAxesFontName','Times')
set(20,'DefaultAxesFontSize',12)
set(20,'DefaultAxesFontWeight','Bold')
set(20,'DefaultTextfontName','Times')
set(20,'DefaultTextfontSize',12)
set(20,'DefaultTextfontWeight','Bold')
set(20,'Position',[0 0 1000 600])

uicontrol('Style', 'text', 'String', 'Choice of y','Position', [250 10 90 20]);
uicontrol('Style', 'text', 'String', 'm','Position', [640 10 15 20]);
uicontrol('Style', 'text', 'String', num2str(y(1),'%+1.2e\n'),'Position', [560 10 80 20]);


  subplot('position',[0.075 0.63 0.33 0.33])
  imagesc(x,z,permute(matxylocalfield(:,1,:),[1 3 2])');
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title('modulus')  
colorbar('vert')

  
  subplot('position',[0.65 0.63 0.33 0.33])
  imagesc(x,z,abs(permute(matxylocalfieldx(:,1,:),[1 3 2])'));
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title(' x')  
colorbar('vert')
  
 subplot('position',[0.075 0.15 0.33 0.33])
  imagesc(x,z,abs(permute(matxylocalfieldy(:,1,:),[1 3 2])'));
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title(' y')  
colorbar('vert')


  subplot('position',[0.65 0.15 0.33 0.33])
  imagesc(x,z,abs(permute(matxylocalfieldz(:,1,:),[1 3 2])'));
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title(' z')  
colorbar('vert')

uicontrol('Style', 'slider', 'Min',1,'Max', ny,...
	  'val',1,'sliderstep',[1/(ny-1) 2/(ny-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotlocalfieldy,x,y,z,matxylocalfield,matxylocalfieldx,matxylocalfieldy,matxylocalfieldz});


case 4


figure(20)
set(20,'DefaultAxesFontName','Times')
set(20,'DefaultAxesFontSize',12)
set(20,'DefaultAxesFontWeight','Bold')
set(20,'DefaultTextfontName','Times')
set(20,'DefaultTextfontSize',12)
set(20,'DefaultTextfontWeight','Bold')
set(20,'Position',[0 0 1000 600])

uicontrol('Style', 'text', 'String', 'Choice of x','Position', [250 10 90 20]);
uicontrol('Style', 'text', 'String', 'm','Position', [640 10 15 20]);
uicontrol('Style', 'text', 'String', num2str(x(1),'%+1.2e\n'),'Position', [560 10 80 20]);

  subplot('position',[0.075 0.63 0.33 0.33])
  imagesc(y,z,permute(matxylocalfield(1,:,:),[2 3 1])');
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title('modulus')  
colorbar('vert')

  
  subplot('position',[0.65 0.63 0.33 0.33])
  imagesc(y,z,abs(permute(matxylocalfieldx(1,:,:),[2 3 1])'));
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title(' x')  
colorbar('vert')
  
 subplot('position',[0.075 0.15 0.33 0.33])
  imagesc(y,z,abs(permute(matxylocalfieldy(1,:,:),[2 3 1])'));
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title('y')  
colorbar('vert')


  subplot('position',[0.65 0.15 0.33 0.33])
  imagesc(y,z,abs(permute(matxylocalfieldz(1,:,:),[2 3 1])'));
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title(' z')  
colorbar('vert')

  

uicontrol('Style', 'slider', 'Min',1,'Max', nx,...
	  'val',1,'sliderstep',[1/(nx-1) 2/(nx-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotlocalfieldx,x,y,z,matxylocalfield,matxylocalfieldx,matxylocalfieldy,matxylocalfieldz});


end;
