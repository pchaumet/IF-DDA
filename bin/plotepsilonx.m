function plotepsilonx(heps,event,x,y,z,matxyepsr,matxyepsi)

% Gets the value of the parameter from the slider.
Param = get(heps,'Value');

% Tranform in integer
i = int8(Param);

% Puts the value of the parameter on the GUI.
uicontrol('Style', 'text', 'String', num2str(x(i),'%+1.2e\n'),...
'Position', [560 15 80 20]);

% Plots the Graph.

figure(2)
set(2,'DefaultAxesFontName','Times')
set(2,'DefaultAxesFontSize',12)
set(2,'DefaultAxesFontWeight','Bold')
set(2,'DefaultTextfontName','Times')
set(2,'DefaultTextfontSize',12)
set(2,'DefaultTextfontWeight','Bold')
set(2,'Position',[0 600 1000 500])

 
  subplot(1,2,1)
  imagesc(y,z,permute(matxyepsr(i,:,:),[2 3 1])');
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title('Re(epsilon)')  
colorbar('vert')

  
  subplot(1,2,2)
  imagesc(y,z,permute(matxyepsi(i,:,:),[2 3 1])');
axis xy

shading interp
xlabel('y')
ylabel('z')
axis image
title('Im(epsilon)')  
colorbar('vert')
  
