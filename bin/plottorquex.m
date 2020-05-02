function plottorquex(hlocal,event,yy,zz,x,matxytorquey,matxytorquez,nprint)

% Gets the value of the parameter from the slider.
Param = get(hlocal,'Value');

% Tranform in integer
i = int8(Param);

% Puts the value of the parameter on the GUI.
uicontrol('Style', 'text', 'String', num2str(x(i)),'Position', [560 15 80 20]);

% Plots the Graph.

figure(300)
set(300,'DefaultAxesFontName','Times')
set(300,'DefaultAxesFontSize',12)
set(300,'DefaultAxesFontWeight','Bold')
set(300,'DefaultTextfontName','Times')
set(300,'DefaultTextfontSize',12)
set(300,'DefaultTextfontWeight','Bold')
set(300,'Position',[0 0 1000 600])
  scale=1;

  quiver(yy(i,:,:),zz(i,:,:),matxytorquey(i,:,:),matxytorquez(i,:,:),scale)
  axis equal
xlabel('y')
ylabel('z')

if (nprint == 1)
print('-f300','torque2d','-depsc')
end
