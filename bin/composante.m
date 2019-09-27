function composante(hObject,event,nx,ny,nz,x,y,z,matxyepsilonr,matxyepsiloni)
  val = get(hObject,'Value');

switch val

case 2

figure(2)
set(2,'DefaultAxesFontName','Times')
set(2,'DefaultAxesFontSize',12)
set(2,'DefaultAxesFontWeight','Bold')
set(2,'DefaultTextfontName','Times')
set(2,'DefaultTextfontSize',12)
set(2,'DefaultTextfontWeight','Bold')
set(2,'Position',[0 600 1000 500])

uicontrol('Style', 'text', 'String', num2str(z(1)),...
'Position', [560 15 60 20]);


  subplot(1,2,1)
  imagesc(x,y,matxyepsilonr(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Re(epsilon)')  
colorbar('vert')

  
  subplot(1,2,2)
  imagesc(x,y,matxyepsiloni(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Im(epsilon)')  
colorbar('vert')
  

uicontrol('Style', 'text', 'String', 'Choice of z',...
'Position', [250 15 90 20]);

uicontrol('Style', 'slider', 'Min',1,'Max', nz,...
	  'val',1,'sliderstep',[1/(nz-1) 2/(nz-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotepsilonz,x,y,z,matxyepsilonr,matxyepsiloni});

case 3


% Plots the Graph.

figure(2)
set(2,'DefaultAxesFontName','Times')
set(2,'DefaultAxesFontSize',12)
set(2,'DefaultAxesFontWeight','Bold')
set(2,'DefaultTextfontName','Times')
set(2,'DefaultTextfontSize',12)
set(2,'DefaultTextfontWeight','Bold')
set(2,'Position',[0 600 1000 500])

uicontrol('Style', 'text', 'String', num2str(y(1)),...
'Position', [560 15 60 20]);


  subplot(1,2,1)
  imagesc(x,z,permute(matxyepsilonr(:,1,:),[1 3 2])');
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title('Re(epsilon)')  
colorbar('vert')

  
  subplot(1,2,2)
  imagesc(x,z,permute(matxyepsiloni(:,1,:),[1 3 2])');
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title('Im(epsilon)')  
colorbar('vert')
  


uicontrol('Style', 'text', 'String', 'Choice of y',...
'Position', [250 15 90 20]);

uicontrol('Style', 'slider', 'Min',1,'Max', ny,...
	  'val',1,'sliderstep',[1/(ny-1) 2/(ny-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotepsilony,x,y,z,matxyepsilonr,matxyepsiloni});
case 4

figure(2)
set(2,'DefaultAxesFontName','Times')
set(2,'DefaultAxesFontSize',12)
set(2,'DefaultAxesFontWeight','Bold')
set(2,'DefaultTextfontName','Times')
set(2,'DefaultTextfontSize',12)
set(2,'DefaultTextfontWeight','Bold')
set(2,'Position',[0 600 1000 500])

uicontrol('Style', 'text', 'String', num2str(x(1)),...
'Position', [560 15 60 20]);

  subplot(1,2,1)
  imagesc(y,z,permute(matxyepsilonr(1,:,:),[2 3 1])');
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title('Re(epsilon)')  
colorbar('vert')

  
  subplot(1,2,2)
  imagesc(y,z,permute(matxyepsiloni(1,:,:),[2 3 1])');
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title('Im(epsilon)')  
colorbar('vert')
  

uicontrol('Style', 'text', 'String', 'Choice of x',...
'Position', [250 15 90 20]);

uicontrol('Style', 'slider', 'Min',1,'Max', nx,...
	  'val',1,'sliderstep',[1/(nx-1) 2/(nx-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotepsilonx,x,y,z,matxyepsilonr,matxyepsiloni});

end;
