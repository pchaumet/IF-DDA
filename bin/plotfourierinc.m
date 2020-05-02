function plotfourierinc(hlocal,event,numaper,kxfourier,fourierm,fourierxc,fourieryc,fourierzc,nprint)

val = get(hlocal,'Value');

switch val

case 1

figure(450)

set(450,'DefaultAxesFontName','Times')
set(450,'DefaultAxesFontSize',12)
set(450,'DefaultAxesFontWeight','Bold')
set(450,'DefaultTextfontName','Times')
set(450,'DefaultTextfontSize',12)
set(450,'DefaultTextfontWeight','Bold')
set(450,'Position',[0 0 1000 600])

 subplot('position',[0.1 0.1 0.8 0.8])

imagesc(kxfourier,kxfourier,fourierm.^2')
axis xy  
caxis([ 0 max(max(fourierm.^2))])
axis equal
axis image
hold on
rectangle('Position',[-numaper -numaper 2*numaper 2*numaper],'Curvature',[1 1],'linewidth',2,'edgecolor','red')


shading interp
axis equal
axis image
colorbar
  
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)
%title('Fourier : modulus','Interpreter','latex','Fontsize',18)


case 2

figure(450)

set(450,'DefaultAxesFontName','Times')
set(450,'DefaultAxesFontSize',12)
set(450,'DefaultAxesFontWeight','Bold')
set(450,'DefaultTextfontName','Times')
set(450,'DefaultTextfontSize',12)
set(450,'DefaultTextfontWeight','Bold')
set(450,'Position',[0 0 1000 600])

 subplot('position',[0.1 0.1 0.8 0.8])

imagesc(kxfourier,kxfourier,fourierm')
axis xy  
caxis([ 0 max(max(fourierm))])
axis equal
axis image
hold on
rectangle('Position',[-numaper -numaper 2*numaper 2*numaper],'Curvature',[1 1],'linewidth',2,'edgecolor','red')


shading interp
axis equal
axis image
colorbar
  
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)
%title('Fourier : modulus','Interpreter','latex','Fontsize',18)

case 3

figure(450)

set(450,'DefaultAxesFontName','Times')
set(450,'DefaultAxesFontSize',12)
set(450,'DefaultAxesFontWeight','Bold')
set(450,'DefaultTextfontName','Times')
set(450,'DefaultTextfontSize',12)
set(450,'DefaultTextfontWeight','Bold')
set(450,'Position',[0 0 1000 600])


%suptitle('Fourier : $x$ component','Interpreter','latex','Fontsize',18)
  
subplot(1,2,1)

  imagesc(kxfourier,kxfourier,abs(fourierxc'))
axis xy

caxis([ 0 max(max(abs(fourierxc')))])
hold on
rectangle('Position',[-numaper -numaper 2*numaper 2*numaper],'Curvature',[1 1],'linewidth',2,'edgecolor','red')

shading interp
axis equal
axis image
colorbar
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)
title('Modulus','Interpreter','latex','Fontsize',18)

subplot(1,2,2)

imagesc(kxfourier,kxfourier,angle(fourierxc)')
axis xy
caxis([-pi pi])


hold on
rectangle('Position',[-numaper -numaper 2*numaper 2*numaper],'Curvature',[1 1],'linewidth',2,'edgecolor','red')

shading interp
axis equal
axis image
colorbar
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)
title('Phase Angle','Interpreter','latex','Fontsize',18)

case 4

figure(450)

set(450,'DefaultAxesFontName','Times')
set(450,'DefaultAxesFontSize',12)
set(450,'DefaultAxesFontWeight','Bold')
set(450,'DefaultTextfontName','Times')
set(450,'DefaultTextfontSize',12)
set(450,'DefaultTextfontWeight','Bold')
set(450,'Position',[0 0 1000 600])


%suptitle('Fourier : $y$ component','Interpreter','latex','Fontsize',18)
  
subplot(1,2,1)

  imagesc(kxfourier,kxfourier,abs(fourieryc'))
axis xy
caxis([ 0 max(max(abs(fourieryc')))])
hold on
rectangle('Position',[-numaper -numaper 2*numaper 2*numaper],'Curvature',[1 1],'linewidth',2,'edgecolor','red')

  
shading interp
axis equal
axis image
colorbar
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)
title('Modulus','Interpreter','latex','Fontsize',18)

subplot(1,2,2)

imagesc(kxfourier,kxfourier,angle(fourieryc)')
axis xy
caxis([-pi pi])

hold on
rectangle('Position',[-numaper -numaper 2*numaper 2*numaper],'Curvature',[1 1],'linewidth',2,'edgecolor','red')

shading interp
axis equal
axis image
colorbar
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)
  
title('Phase Angle','Interpreter','latex','Fontsize',18)

case 5

figure(450)

set(450,'DefaultAxesFontName','Times')
set(450,'DefaultAxesFontSize',12)
set(450,'DefaultAxesFontWeight','Bold')
set(450,'DefaultTextfontName','Times')
set(450,'DefaultTextfontSize',12)
set(450,'DefaultTextfontWeight','Bold')
set(450,'Position',[0 0 1000 600])


%suptitle('Fourier : $z$ component','Interpreter','latex','Fontsize',18)
  
subplot(1,2,1)

  imagesc(kxfourier,kxfourier,abs(fourierzc'))
axis xy
 caxis([ 0 max(max(abs(fourierzc')))])
hold on
rectangle('Position',[-numaper -numaper 2*numaper 2*numaper],'Curvature',[1 1],'linewidth',2,'edgecolor','red')

shading interp
axis equal
axis image
colorbar
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)

title('Modulus','Interpreter','latex','Fontsize',18)

subplot(1,2,2)

imagesc(kxfourier,kxfourier,angle(fourierzc)')
axis xy
caxis([-pi pi])
	       
hold on
rectangle('Position',[-numaper -numaper 2*numaper 2*numaper],'Curvature',[1 1],'linewidth',2,'edgecolor','red')
shading interp
axis equal
axis image
colorbar
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)

title('Phase Angle','Interpreter','latex','Fontsize',18)


end;

if (nprint == 1)
print('-f450','fourierinc','-depsc')
end
