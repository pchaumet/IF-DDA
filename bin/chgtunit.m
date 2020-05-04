function [x,y,z,unitmes] = chgtunit(x,y,z)

xmax=max(x);
ymax=max(y);
zmax=max(z);

xyzmax=max(xmax,ymax);
xyzmax=max(xyzmax,zmax);

unitmes=string('[m]');

if (xyzmax <1) 
x=x*1000;
y=z*1000;
z=z*1000;
xyzmax=xyzmax*1000;
unitmes=string('[mm]');
end;

if (xyzmax <1) 
x=x*1000;
y=z*1000;
z=z*1000;
xyzmax=xyzmax*1000;
unitmes=string('[µm]');
end;

if (xyzmax <1) 
x=x*1000;
y=z*1000;
z=z*1000;
unitmes=string('[nm]');
end;

end
