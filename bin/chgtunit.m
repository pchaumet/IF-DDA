function [x,unitmes] = chgtunit(x)

xmax=max(x);

unitmes=string('[m]');

if (xmax <1) 
x=x*1000;
xmax=xmax*1000;
unitmes=string('[mm]');
end;

if (xmax <1) 
x=x*1000;
xmax=xmax*1000;
unitmes=string('[µm]');
end;

if (xmax <1) 
x=x*1000;
unitmes=string('[nm]');
end;

end
