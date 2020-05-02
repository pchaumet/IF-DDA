function plotforcez(hlocal,event,xx,yy,z,matxyforcex,matxyforcey,nprint)

% Gets the value of the parameter from the slider.
Param = get(hlocal,'Value');

% Tranform in integer
k = int8(Param);

% Puts the value of the parameter on the GUI.
uicontrol('Style', 'text', 'String', num2str(z(k),'%+1.2e\n'),...
'Position', [560 15 80 20]);

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

quiver(xx(:,:,k),yy(:,:,k),matxyforcex(:,:,k),matxyforcey(:,:,k),scale)
  
xlabel('x')
ylabel('y')


if (nprint == 1)
print('-f200','force2d','-depsc')
end
