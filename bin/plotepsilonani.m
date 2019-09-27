function plotepsilonani(hObject,event,nx,ny,nz,x,y,z,matepsilonrxx,matepsilonixx,matepsilonrxy,matepsilonixy,matepsilonrxz,matepsilonixz,matepsilonryx,matepsiloniyx,matepsilonryy,matepsiloniyy,matepsilonryz,matepsiloniyz,matepsilonrzx,matepsilonizx,matepsilonrzy,matepsilonizy,matepsilonrzz,matepsilonizz)

val = get(hObject,'Value')

switch val


case 1

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

maxepsr=max(max(max(matepsilonrxx)));
minepsr=min(min(min(matepsilonrxx)));
maxepsi=max(max(max(matepsilonixx)));
minepsi=min(min(min(matepsilonixx)));
  
  subplot(1,2,1)
  imagesc(x,y,matepsilonrxx(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Re(epsilon)')  
colorbar('vert')
if maxepsr == minepsr; caxis([maxepsr-1 maxepsr+1])
;else; caxis([minepsr maxepsr]);end;
  
  subplot(1,2,2)
  imagesc(x,y,matepsilonixx(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Im(epsilon)')  
colorbar('vert')
if maxepsi == minepsi; caxis([maxepsi-1 maxepsi+1])
;else; caxis([minepsi maxepsi]);end;

uicontrol('Style', 'text', 'String', 'Choice of z',...
'Position', [250 15 90 20]);

uicontrol('Style', 'slider', 'Min',1,'Max', nz,...
	  'val',1,'sliderstep',[1/(nz-1) 2/(nz-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotepsilonz,x,y,z,matepsilonrxx,matepsilonixx});

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

maxepsr=max(max(max(matepsilonrxy)));
minepsr=min(min(min(matepsilonrxy)));
maxepsi=max(max(max(matepsilonixy)));
minepsi=min(min(min(matepsilonixy)));
  
  subplot(1,2,1)
  imagesc(x,y,matepsilonrxy(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Re(epsilon)')  
colorbar('vert')
if maxepsr == minepsr; caxis([maxepsr-1 maxepsr+1])
;else; caxis([minepsr maxepsr]);end;
  
  subplot(1,2,2)
  imagesc(x,y,matepsilonixy(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Im(epsilon)')  
colorbar('vert')
if maxepsi == minepsi; caxis([maxepsi-1 maxepsi+1])
;else; caxis([minepsi maxepsi]);end;

uicontrol('Style', 'text', 'String', 'Choice of z',...
'Position', [250 15 90 20]);

uicontrol('Style', 'slider', 'Min',1,'Max', nz,...
	  'val',1,'sliderstep',[1/(nz-1) 2/(nz-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotepsilonz,x,y,z,matepsilonrxy,matepsilonixy});

case 3

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

maxepsr=max(max(max(matepsilonrxz)));
minepsr=min(min(min(matepsilonrxz)));
maxepsi=max(max(max(matepsilonixz)));
minepsi=min(min(min(matepsilonixz)));
  
  subplot(1,2,1)
  imagesc(x,y,matepsilonrxz(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Re(epsilon)')  
colorbar('vert')
if maxepsr == minepsr; caxis([maxepsr-1 maxepsr+1])
;else; caxis([minepsr maxepsr]);end;
  
  subplot(1,2,2)
  imagesc(x,y,matepsilonixz(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Im(epsilon)')  
colorbar('vert')
if maxepsi == minepsi; caxis([maxepsi-1 maxepsi+1])
;else; caxis([minepsi maxepsi]);end;

uicontrol('Style', 'text', 'String', 'Choice of z',...
'Position', [250 15 90 20]);

uicontrol('Style', 'slider', 'Min',1,'Max', nz,...
	  'val',1,'sliderstep',[1/(nz-1) 2/(nz-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotepsilonz,x,y,z,matepsilonrxz,matepsilonixz});



case 4

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

maxepsr=max(max(max(matepsilonryx)));
minepsr=min(min(min(matepsilonryx)));
maxepsi=max(max(max(matepsiloniyx)));
minepsi=min(min(min(matepsiloniyx)));
  
  subplot(1,2,1)
  imagesc(x,y,matepsilonryx(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Re(epsilon)')  
colorbar('vert')
if maxepsr == minepsr; caxis([maxepsr-1 maxepsr+1])
;else; caxis([minepsr maxepsr]);end;
  
  subplot(1,2,2)
  imagesc(x,y,matepsiloniyx(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Im(epsilon)')  
colorbar('vert')
if maxepsi == minepsi; caxis([maxepsi-1 maxepsi+1])
;else; caxis([minepsi maxepsi]);end;

uicontrol('Style', 'text', 'String', 'Choice of z',...
'Position', [250 15 90 20]);

uicontrol('Style', 'slider', 'Min',1,'Max', nz,...
	  'val',1,'sliderstep',[1/(nz-1) 2/(nz-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotepsilonz,x,y,z,matepsilonryx,matepsiloniyx});

case 5

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

maxepsr=max(max(max(matepsilonryy)));
minepsr=min(min(min(matepsilonryy)));
maxepsi=max(max(max(matepsiloniyy)));
minepsi=min(min(min(matepsiloniyy)));
  
  subplot(1,2,1)
  imagesc(x,y,matepsilonryy(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Re(epsilon)')  
colorbar('vert')
if maxepsr == minepsr; caxis([maxepsr-1 maxepsr+1])
;else; caxis([minepsr maxepsr]);end;
  
  subplot(1,2,2)
  imagesc(x,y,matepsiloniyy(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Im(epsilon)')  
colorbar('vert')
if maxepsi == minepsi; caxis([maxepsi-1 maxepsi+1])
;else; caxis([minepsi maxepsi]);end;

uicontrol('Style', 'text', 'String', 'Choice of z',...
'Position', [250 15 90 20]);

uicontrol('Style', 'slider', 'Min',1,'Max', nz,...
	  'val',1,'sliderstep',[1/(nz-1) 2/(nz-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotepsilonz,x,y,z,matepsilonryy,matepsiloniyy});

case 6

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

maxepsr=max(max(max(matepsilonryz)));
minepsr=min(min(min(matepsilonryz)));
maxepsi=max(max(max(matepsiloniyz)));
minepsi=min(min(min(matepsiloniyz)));
  
  subplot(1,2,1)
  imagesc(x,y,matepsilonryz(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Re(epsilon)')  
colorbar('vert')
if maxepsr == minepsr; caxis([maxepsr-1 maxepsr+1])
;else; caxis([minepsr maxepsr]);end;
  
  subplot(1,2,2)
  imagesc(x,y,matepsiloniyz(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Im(epsilon)')  
colorbar('vert')
if maxepsi == minepsi; caxis([maxepsi-1 maxepsi+1])
;else; caxis([minepsi maxepsi]);end;

uicontrol('Style', 'text', 'String', 'Choice of z',...
'Position', [250 15 90 20]);

uicontrol('Style', 'slider', 'Min',1,'Max', nz,...
	  'val',1,'sliderstep',[1/(nz-1) 2/(nz-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotepsilonz,x,y,z,matepsilonryz,matepsiloniyz});



case 7

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

maxepsr=max(max(max(matepsilonrzx)));
minepsr=min(min(min(matepsilonrzx)));
maxepsi=max(max(max(matepsilonizx)));
minepsi=min(min(min(matepsilonizx)));
  
  subplot(1,2,1)
  imagesc(x,y,matepsilonrzx(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Re(epsilon)')  
colorbar('vert')
if maxepsr == minepsr; caxis([maxepsr-1 maxepsr+1])
;else; caxis([minepsr maxepsr]);end;
  
  subplot(1,2,2)
  imagesc(x,y,matepsilonizx(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Im(epsilon)')  
colorbar('vert')
if maxepsi == minepsi; caxis([maxepsi-1 maxepsi+1])
;else; caxis([minepsi maxepsi]);end;

uicontrol('Style', 'text', 'String', 'Choice of z',...
'Position', [250 15 90 20]);

uicontrol('Style', 'slider', 'Min',1,'Max', nz,...
	  'val',1,'sliderstep',[1/(nz-1) 2/(nz-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotepsilonz,x,y,z,matepsilonrzx,matepsilonizx});

case 8

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

maxepsr=max(max(max(matepsilonrzy)));
minepsr=min(min(min(matepsilonrzy)));
maxepsi=max(max(max(matepsilonizy)));
minepsi=min(min(min(matepsilonizy)));
  
  subplot(1,2,1)
  imagesc(x,y,matepsilonrzy(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Re(epsilon)')  
colorbar('vert')
if maxepsr == minepsr; caxis([maxepsr-1 maxepsr+1])
;else; caxis([minepsr maxepsr]);end;
  
  subplot(1,2,2)
  imagesc(x,y,matepsilonizy(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Im(epsilon)')  
colorbar('vert')
if maxepsi == minepsi; caxis([maxepsi-1 maxepsi+1])
;else; caxis([minepsi maxepsi]);end;

uicontrol('Style', 'text', 'String', 'Choice of z',...
'Position', [250 15 90 20]);

uicontrol('Style', 'slider', 'Min',1,'Max', nz,...
	  'val',1,'sliderstep',[1/(nz-1) 2/(nz-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotepsilonz,x,y,z,matepsilonrzy,matepsilonizy});

case 9

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

maxepsr=max(max(max(matepsilonrzz)));
minepsr=min(min(min(matepsilonrzz)));
maxepsi=max(max(max(matepsilonizz)));
minepsi=min(min(min(matepsilonizz)));
  
  subplot(1,2,1)
  imagesc(x,y,matepsilonrzz(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Re(epsilon)')  
colorbar('vert')
if maxepsr == minepsr; caxis([maxepsr-1 maxepsr+1])
;else; caxis([minepsr maxepsr]);end;
  
  subplot(1,2,2)
  imagesc(x,y,matepsilonizz(:,:,1)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('Im(epsilon)')  
colorbar('vert')
if maxepsi == minepsi; caxis([maxepsi-1 maxepsi+1])
;else; caxis([minepsi maxepsi]);end;

uicontrol('Style', 'text', 'String', 'Choice of z',...
'Position', [250 15 90 20]);

uicontrol('Style', 'slider', 'Min',1,'Max', nz,...
	  'val',1,'sliderstep',[1/(nz-1) 2/(nz-1)],...
	  'Position', [350 10 200 30],'Callback', {@plotepsilonz,x,y,z,matepsilonrzz,matepsilonizz});


end;
