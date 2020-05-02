function plotforcey(hlocal,event,xx,zz,y,matxyforcex,matxyforcez,nprint)

% Gets the value of the parameter from the slider.
Param = get(hlocal,'Value');

% Tranform in integer
j = int8(Param);

% Puts the value of the parameter on the GUI.
uicontrol('Style', 'text', 'String', num2str(y(j),'%+1.2e\n'),'Position', [560 15 80 20]);

% Plots the Graph.

figure(200)
set(200,'DefaultAxesFontName','Times')
set(200,'DefaultAxesFontSize',12)
set(200,'DefaultAxesFontWeight','Bold')
set(200,'DefaultTextfontName','Times')
set(200,'DefaultTextfontSize',12)
set(200,'DefaultTextfontWeight','Bold')
set(200,'Position',[0 0 1000 600])
scale=1

 quiver(xx(:,j,:),zz(:,j,:),matxyforcex(:,j,:),matxyforcez(:,j,:),scale)
  
xlabel('x')
ylabel('z')


if (nprint == 1)
print('-f200','force2d','-depsc')
end
