%% Simple GUI
function SimpleGUIWithoutButton()
close all

load x.mat -ascii
load y.mat -ascii
load z.mat -ascii

  
% Make a large figure.
figure('position',[0 0 700 500], 'name', 'SimpleGUI', 'NumberTitle', 'off');

% Make subplot to hold plot.
h = subplot('position',[0.1 0.3 0.8 0.6]);

% Just some descriptive text.
uicontrol('Style', 'text', 'String', 'Choice of z',...
'Position', [150 50 90 30]);

% A slider for varying the parameter.
uicontrol('Style', 'slider', 'Min',1,'Max', 10,...
 'val',1,'sliderstep',[1/9 2/9],...
'Position', [250 50 200 30],'Callback', @PlotGUI);

%% Called by SimpleGUI to do the plotting
% hObject is the slider and eventdata is unused.
function PlotGUI(hObject,x,y,z)

% Gets the value of the parameter from the slider.
Param = get(hObject,'Value');
load x.mat -ascii
load y.mat -ascii
load z.mat -ascii
load localfield.mat -ascii
matxy=reshape(localfield,10,10,10);
k = Param;

% Puts the value of the parameter on the GUI.
uicontrol('Style', 'text', 'String', num2str(z(k)),...
'Position', [460 55 60 20]);

% Plots the Graph.



imagesc(x,y,matxy(:,:,k));
shading interp
xlabel('x')
ylabel('y')
axis image
title('Local field')  
colorbar('vert')
