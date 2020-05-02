function plotimage(hlocal,event,ximage,imagem,imagexc,imageyc,imagezc,nprint)

val = get(hlocal,'Value');

switch val


case 1


figure(500)

set(500,'DefaultAxesFontName','Times')
set(500,'DefaultAxesFontSize',12)
set(500,'DefaultAxesFontWeight','Bold')
set(500,'DefaultTextfontName','Times')
set(500,'DefaultTextfontSize',12)
set(500,'DefaultTextfontWeight','Bold')
set(500,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.75 0.75])

imagesc(ximage,ximage,imagem.^2')
axis xy
caxis([ 0 max(max(imagem.^2))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)



case 2


figure(500)

set(500,'DefaultAxesFontName','Times')
set(500,'DefaultAxesFontSize',12)
set(500,'DefaultAxesFontWeight','Bold')
set(500,'DefaultTextfontName','Times')
set(500,'DefaultTextfontSize',12)
set(500,'DefaultTextfontWeight','Bold')
set(500,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.75 0.75])

imagesc(ximage,ximage,imagem')
axis xy
caxis([ 0 max(max(imagem))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)


case 3

figure(500)

set(500,'DefaultAxesFontName','Times')
set(500,'DefaultAxesFontSize',12)
set(500,'DefaultAxesFontWeight','Bold')
set(500,'DefaultTextfontName','Times')
set(500,'DefaultTextfontSize',12)
set(500,'DefaultTextfontWeight','Bold')
set(500,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.75 0.75])


  imagesc(ximage,ximage,(imagexc'))
axis xy

caxis([ 0 max(max(abs(imagexc')))])

shading interp
axis equal
axis image
colorbar
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)
%title('Modulus','Interpreter','latex','Fontsize',18)

case 4

figure(500)

set(500,'DefaultAxesFontName','Times')
set(500,'DefaultAxesFontSize',12)
set(500,'DefaultAxesFontWeight','Bold')
set(500,'DefaultTextfontName','Times')
set(500,'DefaultTextfontSize',12)
set(500,'DefaultTextfontWeight','Bold')
set(500,'Position',[0 0 1000 600])


subplot('position',[0.1 0.1 0.75 0.75])

  imagesc(ximage,ximage,(imageyc'))
axis xy
caxis([ 0 max(max(abs(imageyc')))])
  
shading interp
axis equal
axis image
colorbar
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)
%title('Modulus','Interpreter','latex','Fontsize',18)


case 5

figure(500)

set(500,'DefaultAxesFontName','Times')
set(500,'DefaultAxesFontSize',12)
set(500,'DefaultAxesFontWeight','Bold')
set(500,'DefaultTextfontName','Times')
set(500,'DefaultTextfontSize',12)
set(500,'DefaultTextfontWeight','Bold')
set(500,'Position',[0 0 1000 600])


%suptitle('Image : $z$ component','Interpreter','latex','Fontsize',18)
  

subplot('position',[0.1 0.1 0.75 0.75])

  imagesc(ximage,ximage,(imagezc'))

  axis xy
caxis([ 0 max(max((imagezc')))])

shading interp
axis equal
axis image
colorbar
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)

%title('Modulus','Interpreter','latex','Fontsize',18)


end;

if (nprint == 1)
print('-f500','imagewf','-depsc')
end
