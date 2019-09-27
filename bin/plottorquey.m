function plottorquey(hlocal,event,xx,zz,y,matxytorquex,matxytorquez)

% Gets the value of the parameter from the slider.
Param = get(hlocal,'Value');

% Tranform in integer
j = int8(Param);

% Puts the value of the parameter on the GUI.
uicontrol('Style', 'text', 'String', num2str(y(j)),'Position', [560 15 80 20]);

% Plots the Graph.

figure(300)
set(300,'DefaultAxesFontName','Times')
set(300,'DefaultAxesFontSize',12)
set(300,'DefaultAxesFontWeight','Bold')
set(300,'DefaultTextfontName','Times')
set(300,'DefaultTextfontSize',12)
set(300,'DefaultTextfontWeight','Bold')
set(300,'Position',[0 0 1000 600])
scale=1

 quiver(xx(:,j,:),zz(:,j,:),matxytorquex(:,j,:),matxytorquez(:,j,:),scale)
  
xlabel('x')
ylabel('z')

