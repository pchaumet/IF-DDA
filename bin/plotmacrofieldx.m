function plotmacrofieldx(hmacro,event,x,y,z,matxymacrofield,matxymacrofieldx,matxymacrofieldy,matxymacrofieldz,nprint)

% Gets the value of the parameter from the slider.
Param = get(hmacro,'Value');

% Tranform in integer
i = int8(Param);

% Puts the value of the parameter on the GUI.
uicontrol('Style', 'text', 'String', num2str(x(i),'%+1.2e\n'),...
'Position', [560 10 80 20]);

% Plots the Graph.

figure(30)
set(30,'DefaultAxesFontName','Times')
set(30,'DefaultAxesFontSize',12)
set(30,'DefaultAxesFontWeight','Bold')
set(30,'DefaultTextfontName','Times')
set(30,'DefaultTextfontSize',12)
set(30,'DefaultTextfontWeight','Bold')
set(30,'Position',[0 0 1000 600])

  subplot('position',[0.075 0.63 0.33 0.33])
  imagesc(y,z,permute(matxymacrofield(i,:,:),[2 3 1])');
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title('modulus')  
colorbar('vert')

  
  subplot('position',[0.65 0.63 0.33 0.33])
  imagesc(y,z,abs(permute(matxymacrofieldx(i,:,:),[2 3 1])'));
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title(' x')  
colorbar('vert')
  
 subplot('position',[0.075 0.15 0.33 0.33])
  imagesc(y,z,abs(permute(matxymacrofieldy(i,:,:),[2 3 1])'));
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title('y')  
colorbar('vert')

%text(-0.4,0.8,'Macro field','units','Normalized','rotation',90)


  subplot('position',[0.65 0.15 0.33 0.33])
  imagesc(x,y,abs(permute(matxymacrofieldz(i,:,:),[2 3 1])'));
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title(' z')  
colorbar('vert')

if (nprint == 1)
print('-f30','macroscopic','-depsc')
end
