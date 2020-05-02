function plotincifieldy(hinci,event,x,y,z,matxyincifield,matxyincifieldx,matxyincifieldy,matxyincifieldz,nprint)

% Gets the value of the parameter from the slider.
Param = get(hinci,'Value');

% Tranform in integer
j = int8(Param);

% Puts the value of the parameter on the GUI.
uicontrol('Style', 'text', 'String', num2str(y(j),'%+1.2e\n'),...
'Position', [560 10 80 20]);

% Plots the Graph.

figure(10)
set(10,'DefaultAxesFontName','Times')
set(10,'DefaultAxesFontSize',12)
set(10,'DefaultAxesFontWeight','Bold')
set(10,'DefaultTextfontName','Times')
set(10,'DefaultTextfontSize',12)
set(10,'DefaultTextfontWeight','Bold')
set(10,'Position',[0 0 1000 600])

  subplot('position',[0.075 0.63 0.33 0.33])
  imagesc(x,z,permute(matxyincifield(:,j,:),[1 3 2])');
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title('modulus')  
colorbar('vert')

  
  subplot('position',[0.65 0.63 0.33 0.33])
  imagesc(x,z,abs(permute(matxyincifieldx(:,j,:),[1 3 2])'));
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title(' x')  
colorbar('vert')
  
 subplot('position',[0.075 0.15 0.33 0.33])
  imagesc(x,z,abs(permute(matxyincifieldy(:,j,:),[1 3 2])'));
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title(' y')  
colorbar('vert')


  subplot('position',[0.65 0.15 0.33 0.33])
  imagesc(x,z,abs(permute(matxyincifieldz(:,j,:),[1 3 2])'));
axis xy

shading interp
xlabel('x')
ylabel('z')
axis image
title(' z')  
colorbar('vert')

if (nprint == 1)
print('-f10','incident','-depsc')
end

  
