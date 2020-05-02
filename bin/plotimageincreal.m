function plotimageinc(hlocal,event,ximage,imagem,imagexc,imageyc,imagezc,nprint)

val = get(hlocal,'Value');

switch val

case 1


figure(550)

set(550,'DefaultAxesFontName','Times')
set(550,'DefaultAxesFontSize',12)
set(550,'DefaultAxesFontWeight','Bold')
set(550,'DefaultTextfontName','Times')
set(550,'DefaultTextfontSize',12)
set(550,'DefaultTextfontWeight','Bold')
set(550,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.75 0.75])

imagesc(ximage,ximage,imagem.^2')
axis xy
caxis([min(min(imagem.^2)) max(max(imagem.^2))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)


case 2


figure(550)

set(550,'DefaultAxesFontName','Times')
set(550,'DefaultAxesFontSize',12)
set(550,'DefaultAxesFontWeight','Bold')
set(550,'DefaultTextfontName','Times')
set(550,'DefaultTextfontSize',12)
set(550,'DefaultTextfontWeight','Bold')
set(550,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.75 0.75])

imagesc(ximage,ximage,imagem')
axis xy
caxis([min(min(imagem)) max(max(imagem))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)


case 3

figure(550)

set(550,'DefaultAxesFontName','Times')
set(550,'DefaultAxesFontSize',12)
set(550,'DefaultAxesFontWeight','Bold')
set(550,'DefaultTextfontName','Times')
set(550,'DefaultTextfontSize',12)
set(550,'DefaultTextfontWeight','Bold')
set(550,'Position',[0 0 1000 600])


subplot('position',[0.1 0.1 0.75 0.75])


  imagesc(ximage,ximage,(imagexc'))
axis xy

caxis([min(min((imagexc')))  max(max((imagexc')))])

shading interp
axis equal
axis image
colorbar
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)

case 4

figure(550)

set(550,'DefaultAxesFontName','Times')
set(550,'DefaultAxesFontSize',12)
set(550,'DefaultAxesFontWeight','Bold')
set(550,'DefaultTextfontName','Times')
set(550,'DefaultTextfontSize',12)
set(550,'DefaultTextfontWeight','Bold')
set(550,'Position',[0 0 1000 600])


subplot('position',[0.1 0.1 0.75 0.75])


  imagesc(ximage,ximage,(imageyc'))
axis xy
caxis([min(min((imageyc'))) max(max((imageyc')))])
  
shading interp
axis equal
axis image
colorbar
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)

case 5

figure(550)

set(550,'DefaultAxesFontName','Times')
set(550,'DefaultAxesFontSize',12)
set(550,'DefaultAxesFontWeight','Bold')
set(550,'DefaultTextfontName','Times')
set(550,'DefaultTextfontSize',12)
set(550,'DefaultTextfontWeight','Bold')
set(550,'Position',[0 0 1000 600])


subplot('position',[0.1 0.1 0.75 0.75])

  imagesc(ximage,ximage,(imagezc'))

  axis xy
caxis([min(min((imagezc'))) max(max((imagezc')))])

shading interp
axis equal
axis image
colorbar
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)


end;

if (nprint == 1)
print('-f550','imageincwf','-depsc')
end
