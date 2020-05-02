function plotlocalfieldz(hlocal,event,x,y,z,matxylocalfield,matxylocalfieldx,matxylocalfieldy,matxylocalfieldz,nprint)

% Gets the value of the parameter from the slider.
Param = get(hlocal,'Value');

% Tranform in integer
k = int8(Param);

% Puts the value of the parameter on the GUI.
uicontrol('Style', 'text', 'String', num2str(z(k),'%1.2e\n'),...
'Position', [560 10 80 20]);

% Plots the Graph.

figure(20)
set(20,'DefaultAxesFontName','Times')
set(20,'DefaultAxesFontSize',12)
set(20,'DefaultAxesFontWeight','Bold')
set(20,'DefaultTextfontName','Times')
set(20,'DefaultTextfontSize',12)
set(20,'DefaultTextfontWeight','Bold')
set(20,'Position',[0 0 1000 600])

%maxlocalz=max(max(max(matxylocalfield)));
%minlocalz=min(min(min(matxylocalfield)));

  subplot('position',[0.075 0.63 0.33 0.33])
  imagesc(x,y,matxylocalfield(:,:,k)');
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title('modulus')  
colorbar('vert')
%caxis([minlocalz maxlocalz])
  
  subplot('position',[0.65 0.63 0.33 0.33])
  imagesc(x,y,abs(matxylocalfieldx(:,:,k)'));
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title(' x')  
colorbar('vert')

 subplot('position',[0.075 0.15 0.33 0.33])
  imagesc(x,y,abs(matxylocalfieldy(:,:,k)'));
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title(' y')  
colorbar('vert')

%text(-0.4,0.8,'Local field','units','Normalized','rotation',90)


  subplot('position',[0.65 0.15 0.33 0.33])
  imagesc(x,y,abs(matxylocalfieldz(:,:,k)'));
axis xy

shading interp
xlabel('x')
ylabel('y')
axis image
title(' z')  
colorbar('vert')

if (nprint == 1)
print('-f20','local','-depsc')
end
