close all
clear all
nexist=exist('inputmatlab.mat','file');

if (nexist == 0);
namefileh5=input('Please give the name of an H5 file : ','s')
nexisth5 = exist(namefileh5,'file');

if (nexisth5 == 0);
disp('Data files do not exist!')
disp('Please use Advanced interface')
disp('and uncheck the option "Do not write file" if necessary')
return;

end;
end;

if (nexist == 0);

nproche=h5read(namefileh5,'/Option/nproche');
nlocal=h5read(namefileh5,'/Option/nlocal');
nmacro=h5read(namefileh5,'/Option/nmacro');
nsection=h5read(namefileh5,'/Option/nsection');
nsectionsca=h5read(namefileh5,'/Option/nsectionsca');
nquickdiffracte=h5read(namefileh5,'/Option/nquickdiffracte');
nforce=h5read(namefileh5,'/Option/nforce');
nforced=h5read(namefileh5,'/Option/nforced');
ntorque=h5read(namefileh5,'/Option/ntorque');
ntorqued=h5read(namefileh5,'/Option/ntorqued');
nlentille=h5read(namefileh5,'/Option/nlentille');
nquicklens=h5read(namefileh5,'/Option/nquicklens');
nphi=h5read(namefileh5,'/Option/nphi');
ntheta=h5read(namefileh5,'/Option/ntheta');
niso=h5read(namefileh5,'/Option/iso');
nfft=h5read(namefileh5,'/Option/nfft2d');
k0=h5read(namefileh5,'/Option/k0');
numaper=h5read(namefileh5,'/Option/numaper');
nprochefft=h5read(namefileh5,'/Option/nprochefft');
nobjet=h5read(namefileh5,'/Option/nobjet');
ntypemic=h5read(namefileh5,'/Option/ntypemic');
nside=h5read(namefileh5,'/Option/nside');
ntypefile=h5read(namefileh5,'/Option/nmat');
numaperinc=h5read(namefileh5,'/Option/numaperinc');

else

load inputmatlab.mat -ascii
nproche=inputmatlab(1);         % Defined box of computation for the near field
nlocal=inputmatlab(2);          % Compute the local field
nmacro=inputmatlab(3);          % Compute the macroscopic field
nsection=inputmatlab(4);        % Compute the cross section
nsectionsca=inputmatlab(5);     % Compute  C_sca, Poynting and g
nquickdiffracte=inputmatlab(6); % Compute  C_sca, Poynting and g with FFT
nforce=inputmatlab(7);          % Compute the optical force
nforced=inputmatlab(8);         % Compute the optical force
ntorque=inputmatlab(9);         % Compute the optical force
ntorqued=inputmatlab(10);       % Compute the optical force
nlentille=inputmatlab(11);      % Compute the object through the microscope
nquicklens=inputmatlab(12);     % Compte the lens with FFT
nphi=inputmatlab(13);           % nphi
ntheta=inputmatlab(14);         % ntheta
niso=inputmatlab(15);           % 0 isotrope, 1 anisotrope
nfft=inputmatlab(16);           % size of the FFT
k0=inputmatlab(17);             % Wavenumber
numaper=inputmatlab(18);        % Numerical aperture
nprochefft=inputmatlab(19);     % =1 if wide field
nobjet=inputmatlab(20);         % =1 do only the objet
ntypemic=inputmatlab(21);       % type de microscope
nside=inputmatlab(22);          % reflexion or transmission
ntypefile=inputmatlab(23);      % mat file or hdf5 file
numaperinc=inputmatlab(24);     % Numerical aperture for condenser

if (ntypefile == 0);
disp('Use mat file')
end;


if (ntypefile == 1);
disp('Data files do not computed')
disp('Please use Advanced interface')
disp('and uncheck the option "Do not write file" if necessary')
return;
end;


if (ntypefile == 2);


nexisth5a = exist('filenameh5','file');
if (nexisth5a == 0);
disp('H5 files do not exist!')
return;
end;

fid=fopen('filenameh5');    % name file of hdf5
namefileh5a=fgetl(fid);
k=0; for i=1:101; if namefileh5a(i) ~= ' '; k = k+1; namefileh5(k)=namefileh5a(i);end;end;
disp('Use HDF5 file')
%namefileh5=input('Please give the name of an H5 file : ','s')
end;

end;

icomp=complex(0,1);

nprint = input('Print figures in eps (yes=1)')


%%%%%%%%%%%%%%%% Begin plot dipole %%%%%%%%%%%%%%%%%%%%%%%%%%

if (ntypefile == 0);
load x.mat -ascii
load y.mat -ascii
load z.mat -ascii

nx=max(size(x));
ny=max(size(y));
nz=max(size(z));

load xc.mat -ascii
load yc.mat -ascii
load zc.mat -ascii
load epsilon.mat -ascii

elseif(ntypefile == 2)

nx=double(h5read(namefileh5,'/Object/nx'));
ny=double(h5read(namefileh5,'/Object/ny'));
nz=double(h5read(namefileh5,'/Object/nz'));
xc=h5read(namefileh5,'/Object/Dipole position x');
yc=h5read(namefileh5,'/Object/Dipole position y');
zc=h5read(namefileh5,'/Object/Dipole position z');
epsilon(:,1)=h5read(namefileh5,'/Object/Epsilon real part');
epsilon(:,2)=h5read(namefileh5,'/Object/Epsilon imaginary part');

xmin=min(xc);
xmax=max(xc);
pasx=(xmax-xmin)/(nx-1);
ymin=min(yc);
ymax=max(yc);
pasy=(ymax-ymin)/(ny-1);
zmin=min(zc);
zmax=max(zc);
pasz=(zmax-zmin)/(nz-1);
x=xmin:pasx:xmax;
y=ymin:pasy:ymax;
z=zmin:pasz:zmax;

end;
  

if (niso == 0);
if (nproche == -1 );

figure(1)
set(1,'DefaultAxesFontName','Times')
set(1,'DefaultAxesFontSize',12)
set(1,'DefaultAxesFontWeight','Bold')
set(1,'DefaultTextfontName','Times')
set(1,'DefaultTextfontSize',12)
set(1,'DefaultTextfontWeight','Bold')
set(1,'Position',[0 600 500 500])


m=0;n=max(size( epsilon)); 
for i=1:n; if epsilon(i,1)~= 1 || epsilon(i,2)~= 0;
m=m+1;epsilonbg(m)=epsilon(i,1);end;end;

  
scatter3(xc,yc,zc,10,epsilonbg)
axis image
colorbar  
xlabel('x')
ylabel('y')
zlabel('z')  
title('Position of the dipoles')

if (nprint == 1)
print('-f1','dipolepos','-depsc')
end

if nobjet == 1; return ;end;

 else

figure(1)
set(1,'DefaultAxesFontName','Times')
set(1,'DefaultAxesFontSize',12)
set(1,'DefaultAxesFontWeight','Bold')
set(1,'DefaultTextfontName','Times')
set(1,'DefaultTextfontSize',12)
set(1,'DefaultTextfontWeight','Bold')
set(1,'Position',[0 600 1200 500])

  subplot(1,2,1)
     
scatter3(xc,yc,zc,10,epsilon(:,1))
axis image
colorbar  
xlabel('x')
ylabel('y')
zlabel('z')  
title('Position of the dipoles')
subplot(1,2,2)

m=0;n=max(size( epsilon)); 
for i=1:n; if epsilon(i,1)~= 1 || epsilon(i,2)~= 0 ; m=m+1;
xcc(m)=xc(i);ycc(m)=yc(i);zcc(m)=zc(i);epsilonbg(m)=epsilon(i,1);
end;end;


  
scatter3(xcc,ycc,zcc,10,epsilonbg)
axis image
colorbar  
xlabel('x')
ylabel('y')
zlabel('z')  
title({'Position of the dipoles' 'without background'})
end;

if (nprint == 1)
print('-f1','dipolepos','-depsc')
end

else

figure(1)
set(1,'DefaultAxesFontName','Times')
set(1,'DefaultAxesFontSize',12)
set(1,'DefaultAxesFontWeight','Bold')
set(1,'DefaultTextfontName','Times')
set(1,'DefaultTextfontSize',12)
set(1,'DefaultTextfontWeight','Bold')
set(1,'Position',[0 600 1200 500])

 
     
  scatter3(xc,yc,zc,10)
  axis image
colorbar  
xlabel('x')
ylabel('y')
zlabel('z')  
title('Position of the dipoles (anisotropic object)')

if (nprint == 1)
print('-f1','dipolepos','-depsc')
end

if nobjet == 1; return ;end;
end;
%%%%%%%%%%%%%%%% End plot dipole %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Begin plot epsilon %%%%%%%%%%%%%%%%%%%%%%%%%%



if niso == 0;


matxyepsilonr=reshape(epsilon(:,1),nx,ny,nz);
matxyepsiloni=reshape(epsilon(:,2),nx,ny,nz);

clear epsilon

figure(2)
set(2,'DefaultAxesFontName','Times')
set(2,'DefaultAxesFontSize',12)
set(2,'DefaultAxesFontWeight','Bold')
set(2,'DefaultTextfontName','Times')
set(2,'DefaultTextfontSize',12)
set(2,'DefaultTextfontWeight','Bold')
set(2,'Position',[0 600 1000 500])

uicontrol('Style','text','Fontsize',16,'Fontweight','bold',...
'Position',[380 440 300 50],'String','Plot relative permittivity');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 400 170 50],...
	  'Callback',{@plotepsilon,nx,ny,nz,x,y,z,matxyepsilonr,matxyepsiloni,nprint});

if (nprint == 1)
print('-f2','epsilon','-depsc')
end

else

for i=1:nz;for j=1:ny;for k=1:nx;kk=9*(k+nx*(j-1)+nx*ny*(i-1)-1);
matepsilonrxx(i,j,k)=epsilon(kk+1,1);
matepsilonixx(i,j,k)=epsilon(kk+1,2);
matepsilonrxy(i,j,k)=epsilon(kk+2,1);
matepsilonixy(i,j,k)=epsilon(kk+2,2);
matepsilonrxz(i,j,k)=epsilon(kk+3,1);
matepsilonixz(i,j,k)=epsilon(kk+3,2);
matepsilonryx(i,j,k)=epsilon(kk+4,1);
matepsiloniyx(i,j,k)=epsilon(kk+4,2);
matepsilonryy(i,j,k)=epsilon(kk+5,1);
matepsiloniyy(i,j,k)=epsilon(kk+5,2);
matepsilonryz(i,j,k)=epsilon(kk+6,1);
matepsiloniyz(i,j,k)=epsilon(kk+6,2);
matepsilonrzx(i,j,k)=epsilon(kk+7,1);
matepsilonizx(i,j,k)=epsilon(kk+7,2);
matepsilonrzy(i,j,k)=epsilon(kk+8,1);
matepsilonizy(i,j,k)=epsilon(kk+8,2);
matepsilonrzz(i,j,k)=epsilon(kk+9,1);
matepsilonizz(i,j,k)=epsilon(kk+9,2);
end;end;end;
clear epsilon

  
figure(2)
set(2,'DefaultAxesFontName','Times')
set(2,'DefaultAxesFontSize',12)
set(2,'DefaultAxesFontWeight','Bold')
set(2,'DefaultTextfontName','Times')
set(2,'DefaultTextfontSize',12)
set(2,'DefaultTextfontWeight','Bold')
set(2,'Position',[0 600 1000 500])

uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[300 440 380 50],'String','Plot tensor of permittivity in (x,y)-plane for        component');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
	  {'xx','xy','xz','yx','yy','yz','zx','zy','zz'},...
'Position',[490 425 40 40],...
	  'Callback',{@plotepsilonani,nx,ny,nz,x,y,z,matepsilonrxx,matepsilonixx,matepsilonrxy,matepsilonixy,matepsilonrxz,matepsilonixz,matepsilonryx,matepsiloniyx,matepsilonryy,matepsiloniyy,matepsilonryz,matepsiloniyz,matepsilonrzx,matepsilonizx,matepsilonrzy,matepsilonizy,matepsilonrzz,matepsilonizz});

if (nprint == 1)
print('-f2','epsilon','-depsc')
end
		   
end;
%%%%%%%%%%%%%%%% End plot epsilon %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Begin incident field %%%%%%%%%%%%%%%%%%%%%%%%%%

if ( nprochefft == 0) 

if (ntypefile == 0);
 load incidentfield.mat -ascii
 load incidentfieldx.mat -ascii
 load incidentfieldy.mat -ascii
 load incidentfieldz.mat -ascii
elseif (ntypefile ==2);
incidentfield=h5read(namefileh5,'/Near Field/Incident field modulus');
incidentfieldx(:,1)=h5read(namefileh5,'/Near Field/Incident field x component real part');
incidentfieldx(:,2)=h5read(namefileh5,'/Near Field/Incident field x component imaginary part');
incidentfieldy(:,1)=h5read(namefileh5,'/Near Field/Incident field y component real part');
incidentfieldy(:,2)=h5read(namefileh5,'/Near Field/Incident field y component imaginary part');
incidentfieldz(:,1)=h5read(namefileh5,'/Near Field/Incident field z component real part');
incidentfieldz(:,2)=h5read(namefileh5,'/Near Field/Incident field z component imaginary part');
end;

matxyincifield=reshape(incidentfield,nx,ny,nz);
matxyincifieldx=reshape(incidentfieldx(:,1)+icomp*incidentfieldx(:,2),nx,ny,nz);
matxyincifieldy=reshape(incidentfieldy(:,1)+icomp*incidentfieldy(:,2),nx,ny,nz);
matxyincifieldz=reshape(incidentfieldz(:,1)+icomp*incidentfieldz(:,2),nx,ny,nz);

clear incidentfieldx
clear incidentfieldy
clear incidentfieldz


figure(10)
set(10,'DefaultAxesFontName','Times')
set(10,'DefaultAxesFontSize',12)
set(10,'DefaultAxesFontWeight','Bold')
set(10,'DefaultTextfontName','Times')
set(10,'DefaultTextfontSize',12)
set(10,'DefaultTextfontWeight','Bold')
set(10,'Position',[0 0 1000 600])




uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[400 570 210 20],'String','Plot incident field');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
	  'Callback',{@plotincifield,nx,ny,nz,x,y,z,matxyincifield,matxyincifieldx,matxyincifieldy,matxyincifieldz,nprint});

if (nprint == 1)
print('-f10','incident','-depsc')
end

 else

   
load xwf.mat -ascii
load ywf.mat -ascii
load zwf.mat -ascii

nxm=max(size(xwf));
nym=max(size(ywf));
nzm=max(size(zwf));


if (ntypefile == 0);
load incidentfieldwf.mat -ascii
load incidentfieldxwf.mat -ascii
load incidentfieldywf.mat -ascii
load incidentfieldzwf.mat -ascii
elseif (ntypefile ==2);
incidentfieldwf=h5read(namefileh5,'/Near Field/Incident field modulus wf');
incidentfieldxwf(:,1)=h5read(namefileh5,'/Near Field/Incident field x component real part wf');
incidentfieldxwf(:,2)=h5read(namefileh5,'/Near Field/Incident field x component imaginary part wf');
incidentfieldywf(:,1)=h5read(namefileh5,'/Near Field/Incident field y component real part wf');
incidentfieldywf(:,2)=h5read(namefileh5,'/Near Field/Incident field y component imaginary part wf');
incidentfieldzwf(:,1)=h5read(namefileh5,'/Near Field/Incident field z component real part wf');
incidentfieldzwf(:,2)=h5read(namefileh5,'/Near Field/Incident field z component imaginary part wf');
end;


matxyincifield=reshape(incidentfieldwf,nxm,nym,nzm);
matxyincifieldx=reshape(incidentfieldxwf(:,1)+icomp*incidentfieldxwf(:,2),nxm,nym,nzm);
matxyincifieldy=reshape(incidentfieldywf(:,1)+icomp*incidentfieldywf(:,2),nxm,nym,nzm);
matxyincifieldz=reshape(incidentfieldzwf(:,1)+icomp*incidentfieldzwf(:,2),nxm,nym,nzm);

clear incidentfieldxwf
clear incidentfieldywf
clear incidentfieldzwf

figure(10)
set(10,'DefaultAxesFontName','Times')
set(10,'DefaultAxesFontSize',12)
set(10,'DefaultAxesFontWeight','Bold')
set(10,'DefaultTextfontName','Times')
set(10,'DefaultTextfontSize',12)
set(10,'DefaultTextfontWeight','Bold')
set(10,'Position',[0 0 1000 600])




uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[400 570 210 20],'String','Plot incident field');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
	  'Callback',{@plotincifield,nx,ny,nz,x,y,z,matxyincifield,matxyincifieldx,matxyincifieldy,matxyincifieldz,nprint});

if (nprint == 1)
print('-f10','incident','-depsc')
end

end;



%%%%%%%%%%%%%%%% End incident field %%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% Plot local field %%%%%%%%%%%%%%%%%%%%%%%%%%

  
if nlocal == 1;

if (nprochefft ==0) 

if (ntypefile == 0);  
 load localfield.mat -ascii
 load localfieldx.mat -ascii
 load localfieldy.mat -ascii
 load localfieldz.mat -ascii
elseif (ntypefile ==2);
localfield=h5read(namefileh5,'/Near Field/Local field modulus');
localfieldx(:,1)=h5read(namefileh5,'/Near Field/Local field x component real part');
localfieldx(:,2)=h5read(namefileh5,'/Near Field/Local field x component imaginary part');
localfieldy(:,1)=h5read(namefileh5,'/Near Field/Local field y component real part');
localfieldy(:,2)=h5read(namefileh5,'/Near Field/Local field y component imaginary part');
localfieldz(:,1)=h5read(namefileh5,'/Near Field/Local field z component real part');
localfieldz(:,2)=h5read(namefileh5,'/Near Field/Local field z component imaginary part');
end;

matxylocalfield=reshape(localfield,nx,ny,nz);
matxylocalfieldx=reshape(localfieldx(:,1)+icomp*localfieldx(:,2),nx,ny,nz);
matxylocalfieldy=reshape(localfieldy(:,1)+icomp*localfieldy(:,2),nx,ny,nz);
matxylocalfieldz=reshape(localfieldz(:,1)+icomp*localfieldz(:,2),nx,ny,nz);

clear localfieldx
clear localfieldy
clear localfieldz

figure(20)
set(20,'DefaultAxesFontName','Times')
set(20,'DefaultAxesFontSize',12)
set(20,'DefaultAxesFontWeight','Bold')
set(20,'DefaultTextfontName','Times')
set(20,'DefaultTextfontSize',12)
set(20,'DefaultTextfontWeight','Bold')
set(20,'Position',[0 0 1000 600])



uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[400 570 200 20],'String','Plot local field');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
	  'Callback',{@plotlocalfield,nx,ny,nz,x,y,z,matxylocalfield,matxylocalfieldx,matxylocalfieldy,matxylocalfieldz,nprint});

if (nprint == 1)
print('-f20','local','-depsc')
end

 else
   
load xwf.mat -ascii
load ywf.mat -ascii
load zwf.mat -ascii

nxm=max(size(xwf));
nym=max(size(ywf));
nzm=max(size(zwf));

if (ntypefile ==0);
load localfieldwf.mat -ascii
load localfieldxwf.mat -ascii
load localfieldywf.mat -ascii
load localfieldzwf.mat -ascii
elseif (ntypefile ==2);
localfieldwf=h5read(namefileh5,'/Near Field/Local field modulus wf');
localfieldxwf(:,1)=h5read(namefileh5,'/Near Field/Local field x component real part wf');
localfieldxwf(:,2)=h5read(namefileh5,'/Near Field/Local field x component imaginary part wf');
localfieldywf(:,1)=h5read(namefileh5,'/Near Field/Local field y component real part wf');
localfieldywf(:,2)=h5read(namefileh5,'/Near Field/Local field y component imaginary part wf');
localfieldzwf(:,1)=h5read(namefileh5,'/Near Field/Local field z component real part wf');
localfieldzwf(:,2)=h5read(namefileh5,'/Near Field/Local field z component imaginary part wf');

end;

  
matxylocalfield=reshape(localfieldwf,nxm,nym,nzm);
matxylocalfieldx=reshape(localfieldxwf(:,1)+icomp*localfieldxwf(:,2),nxm,nym,nzm);
matxylocalfieldy=reshape(localfieldywf(:,1)+icomp*localfieldywf(:,2),nxm,nym,nzm);
matxylocalfieldz=reshape(localfieldzwf(:,1)+icomp*localfieldzwf(:,2),nxm,nym,nzm);

clear localfieldxwf
clear localfieldywf
clear localfieldzwf

figure(20)
set(20,'DefaultAxesFontName','Times')
set(20,'DefaultAxesFontSize',12)
set(20,'DefaultAxesFontWeight','Bold')
set(20,'DefaultTextfontName','Times')
set(20,'DefaultTextfontSize',12)
set(20,'DefaultTextfontWeight','Bold')
set(20,'Position',[0 0 1000 600])


uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[400 570 200 20],'String','Plot local field');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
	  'Callback',{@plotlocalfield,nxm,nym,nzm,xwf,ywf,zwf,matxylocalfield,matxylocalfieldx,matxylocalfieldy,matxylocalfieldz,nprint});

if (nprint == 1)
print('-f20','local','-depsc')
end

end;


end;

%%%%%%%%%%%%%%%% End local field %%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%% Begin macrocopic field %%%%%%%%%%%%%%%%%%%%%%%%%%
if nmacro == 1;


if (nprochefft ==0) 
if (ntypefile == 0);
load macroscopicfield.mat -ascii
load macroscopicfieldx.mat -ascii
load macroscopicfieldy.mat -ascii
load macroscopicfieldz.mat -ascii
elseif (ntypefile ==2);
macroscopicfield=h5read(namefileh5,'/Near Field/Macroscopic field modulus');
macroscopicfieldx(:,1)=h5read(namefileh5,'/Near Field/Macroscopic field x component real part');
macroscopicfieldx(:,2)=h5read(namefileh5,'/Near Field/Macroscopic field x component imaginary part');
macroscopicfieldy(:,1)=h5read(namefileh5,'/Near Field/Macroscopic field y component real part');
macroscopicfieldy(:,2)=h5read(namefileh5,'/Near Field/Macroscopic field y component imaginary part');
macroscopicfieldz(:,1)=h5read(namefileh5,'/Near Field/Macroscopic field z component real part');
macroscopicfieldz(:,2)=h5read(namefileh5,'/Near Field/Macroscopic field z component imaginary part');
end;

matxymacrofield=reshape(macroscopicfield,nx,ny,nz);
matxymacrofieldx=reshape(macroscopicfieldx(:,1)+icomp*macroscopicfieldx(:,2),nx,ny,nz);
matxymacrofieldy=reshape(macroscopicfieldy(:,1)+icomp*macroscopicfieldy(:,2),nx,ny,nz);
matxymacrofieldz=reshape(macroscopicfieldz(:,1)+icomp*macroscopicfieldz(:,2),nx,ny,nz);

clear macroscopicfieldx
clear macroscopicfieldy
clear macroscopicfieldz

figure(30)
set(30,'DefaultAxesFontName','Times')
set(30,'DefaultAxesFontSize',12)
set(30,'DefaultAxesFontWeight','Bold')
set(30,'DefaultTextfontName','Times')
set(30,'DefaultTextfontSize',12)
set(30,'DefaultTextfontWeight','Bold')
set(30,'Position',[0 0 1000 600])



uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[400 550 200 50],'String','Plot macroscopic field')

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
	  'Callback',{@plotmacrofield,nx,ny,nz,x,y,z,matxymacrofield,matxymacrofieldx,matxymacrofieldy,matxymacrofieldz,nprint});

if (nprint == 1)
print('-f30','macroscopic','-depsc')
end


 else
load xwf.mat -ascii
load ywf.mat -ascii
load zwf.mat -ascii

nxm=max(size(xwf));
nym=max(size(ywf));
nzm=max(size(zwf));
if (ntypefile == 0);  
load macroscopicfieldwf.mat -ascii
load macroscopicfieldxwf.mat -ascii
load macroscopicfieldywf.mat -ascii
load macroscopicfieldzwf.mat -ascii
elseif (ntypefile ==2);
macroscopicfieldwf=h5read(namefileh5,'/Near Field/Macroscopic field modulus wf');
macroscopicfieldxwf(:,1)=h5read(namefileh5,'/Near Field/Macroscopic field x component real part wf');
macroscopicfieldxwf(:,2)=h5read(namefileh5,'/Near Field/Macroscopic field x component imaginary part wf');
macroscopicfieldywf(:,1)=h5read(namefileh5,'/Near Field/Macroscopic field y component real part wf');
macroscopicfieldywf(:,2)=h5read(namefileh5,'/Near Field/Macroscopic field y component imaginary part wf');
macroscopicfieldzwf(:,1)=h5read(namefileh5,'/Near Field/Macroscopic field z component real part wf');
macroscopicfieldzwf(:,2)=h5read(namefileh5,'/Near Field/Macroscopic field z component imaginary part wf');
end;

  
matxymacrofield=reshape(macroscopicfieldwf,nxm,nym,nzm);
matxymacrofieldx=reshape(macroscopicfieldxwf(:,1)+icomp*macroscopicfieldxwf(:,2),nxm,nym,nzm);
matxymacrofieldy=reshape(macroscopicfieldywf(:,1)+icomp*macroscopicfieldywf(:,2),nxm,nym,nzm);
matxymacrofieldz=reshape(macroscopicfieldzwf(:,1)+icomp*macroscopicfieldzwf(:,2),nxm,nym,nzm);

clear macrofieldxwf
clear macrofieldywf
clear macrofieldzwf

figure(30)
set(30,'DefaultAxesFontName','Times')
set(30,'DefaultAxesFontSize',12)
set(30,'DefaultAxesFontWeight','Bold')
set(30,'DefaultTextfontName','Times')
set(30,'DefaultTextfontSize',12)
set(30,'DefaultTextfontWeight','Bold')
set(30,'Position',[0 0 1000 600])



uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[400 550 200 50],'String','Plot macroscopic field')

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
	  'Callback',{@plotmacrofield,nxm,nym,nzm,xwf,ywf,zwf,matxymacrofield,matxymacrofieldx,matxymacrofieldy,matxymacrofieldz,nprint});

if (nprint == 1)
print('-f30','macroscopic','-depsc')
end

end;



end;
%%%%%%%%%%%%%%%% End macrocopic field %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Begin Poynting vector %%%%%%%%%%%%%%%%%%%%%%%%%%
if nsectionsca == 1;
if nquickdiffracte == 0;

if (ntypefile ==0)
load poynting.mat -ascii
elseif (ntypefile ==2)
poynting=h5read(namefileh5,'/Far Field/Poynting');
end;

ray=reshape(poynting,nphi,ntheta);
ray(nphi+1,:)=ray(1,:);
  
for i=1:ntheta;for j=1:nphi+1;
theta=pi*(i-1)/(ntheta-1);
phi=2*pi*(j-1)/(nphi);
xpo(j,i)=ray(j,i)*sin(theta)*cos(phi);
ypo(j,i)=ray(j,i)*sin(theta)*sin(phi);
zpo(j,i)=ray(j,i)*cos(theta);
end;
end;


figure(100)
set(100,'DefaultAxesFontName','Times')
set(100,'DefaultAxesFontSize',12)
set(100,'DefaultAxesFontWeight','Bold')
set(100,'DefaultTextfontName','Times')
set(100,'DefaultTextfontSize',12)
set(100,'DefaultTextfontWeight','Bold')
set(100,'Position',[1000 0 600 600])
surf(xpo,ypo,zpo,ray)
shading interp

title('Poynting vector')

xlabel('x')
ylabel('y')
zlabel('z')

if (nprint == 1)
print('-f100','poynting3d','-depsc')
end

else


if (ntypefile ==0)
load poynting.mat -ascii
elseif (ntypefile ==2)
poynting=h5read(namefileh5,'/Far Field/Poynting');
end;

ray=reshape(poynting,nphi,ntheta);
ray(nphi+1,:)=ray(1,:);
  
for i=1:ntheta;for j=1:nphi+1;
theta=pi*(i-1)/(ntheta-1);
phi=2*pi*(j-1)/(nphi);
xpo(j,i)=ray(j,i)*sin(theta)*cos(phi);
ypo(j,i)=ray(j,i)*sin(theta)*sin(phi);
zpo(j,i)=ray(j,i)*cos(theta);
end;
end;


figure(100)
set(100,'DefaultAxesFontName','Times')
set(100,'DefaultAxesFontSize',12)
set(100,'DefaultAxesFontWeight','Bold')
set(100,'DefaultTextfontName','Times')
set(100,'DefaultTextfontSize',12)
set(100,'DefaultTextfontWeight','Bold')
set(100,'Position',[1000 0 600 600])

surf(xpo,ypo,zpo,ray)
shading interp

title('Poynting vector interpolated in 3D')

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)
zlabel('$z$','Interpreter','latex','Fontsize',18)

if (nprint == 1)
print('-f100','poynting3d','-depsc')
end

if (ntypefile ==0)
load poyntingpos.mat -ascii
load poyntingneg.mat -ascii
load kx.mat -ascii
load ky.mat -ascii
elseif (ntypefile ==2)
poyntingpos=h5read(namefileh5,'/Far Field/Poynting positive');
poyntingneg=h5read(namefileh5,'/Far Field/Poynting negative');
kx=h5read(namefileh5,'/Far Field/kx Poynting');
ky=h5read(namefileh5,'/Far Field/ky Poynting');
end;

nxp=max(size(kx));
nyp=max(size(ky));



figure(101)
set(101,'DefaultAxesFontName','Times')
set(101,'DefaultAxesFontSize',12)
set(101,'DefaultAxesFontWeight','Bold')
set(101,'DefaultTextfontName','Times')
set(101,'DefaultTextfontSize',12)
set(101,'DefaultTextfontWeight','Bold')
set(101,'Position',[1000 0 1000 600])

suptitle('Poynting Modulus in $k_x$ and $k_y$ plane','Interpreter','latex','Fontsize',18)

subplot(1,2,1)

imagesc(kx/k0,ky/k0,reshape(poyntingpos,nxp,nyp)')
axis xy

hold on
rectangle('Position',[-1 -1 2 2],'Curvature',[1 1],'linewidth',2,'edgecolor','red')

axis image
axis equal
title('$k_z>0$','Interpreter','latex','Fontsize',18)

xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)
colorbar

subplot(1,2,2)

imagesc(kx/k0,ky/k0,reshape(poyntingneg,nxp,nyp)')
axis xy

hold on
rectangle('Position',[-1 -1 2 2],'Curvature',[1 1],'linewidth',2,'edgecolor','red')

axis image
axis equal
title('$k_z<0$','Interpreter','latex','Fontsize',18)

xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)

colorbar

if (nprint == 1)
print('-f101','poynting2d','-depsc')
end


end;

end;
%%%%%%%%%%%%%%%% End Poynting vector %%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%% Begin optical force %%%%%%%%%%%%%%%%%%%%%%%%%%
if nforced==1;

if (ntypefile ==0)
load forcex.mat -ascii
load forcey.mat -ascii
load forcez.mat -ascii
  elseif (ntypefile ==2)
forcex=h5read(namefileh5,'/Optical Force/Optical force x component');
forcey=h5read(namefileh5,'/Optical Force/Optical force y component');
forcez=h5read(namefileh5,'/Optical Force/Optical force z component');
  end;


matxyforcex=reshape(forcex,nx,ny,nz);
matxyforcey=reshape(forcey,nx,ny,nz);
matxyforcez=reshape(forcez,nx,ny,nz);

clear forcex
clear forcey
clear forcez

for i=1:nz;for j=1:ny;for k=1:nx;
xx(k,j,i)=x(k);yy(k,j,i)=y(j);zz(k,j,i)=z(i);
end;end;end;

figure(200)
set(200,'DefaultAxesFontName','Times')
set(200,'DefaultAxesFontSize',12)
set(200,'DefaultAxesFontWeight','Bold')
set(200,'DefaultTextfontName','Times')
set(200,'DefaultTextfontSize',12)
set(200,'DefaultTextfontWeight','Bold')
set(200,'Position',[0 0 1000 600])



uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[350 570 300 25],'String','Plot Optical force');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
	  'Callback',{@plotforce,nx,ny,nz,x,y,z,xx,yy,zz,matxyforcex,matxyforcey,matxyforcez});

if (nprint == 1)
print('-f200','force2d','-depsc')
end

figure(201)
set(201,'DefaultAxesFontName','Times')
set(201,'DefaultAxesFontSize',12)
set(201,'DefaultAxesFontWeight','Bold')
set(201,'DefaultTextfontName','Times')
set(201,'DefaultTextfontSize',12)
set(201,'DefaultTextfontWeight','Bold')
set(201,'Position',[0 0 1000 600])

quiver3(xx,yy,zz,matxyforcex,matxyforcey,matxyforcez)


xlabel('x')
ylabel('y')
zlabel('z')
title('Density of optical force')

if (nprint == 1)
print('-f201','force3d','-depsc')
end

end;  
%%%%%%%%%%%%%%%% End  optical force %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Begin optical torque %%%%%%%%%%%%%%%%%%%%%%%%%%
if ntorqued==1;

if (ntypefile ==0)
load torquex.mat -ascii
load torquey.mat -ascii
load torquez.mat -ascii
elseif (ntypefile ==2)
torquex=h5read(namefileh5,'/Optical Force/Optical torque x component');
torquey=h5read(namefileh5,'/Optical Force/Optical torque y component');
torquez=h5read(namefileh5,'/Optical Force/Optical torque z component');
end;

matxytorquex=reshape(torquex,nx,ny,nz);
matxytorquey=reshape(torquey,nx,ny,nz);
matxytorquez=reshape(torquez,nx,ny,nz);

clear torquex
clear torquey
clear torquez

for i=1:nz;for j=1:ny;for k=1:nx;
xx(k,j,i)=x(k);yy(k,j,i)=y(j);zz(k,j,i)=z(i);
end;end;end;

figure(300)
set(300,'DefaultAxesFontName','Times')
set(300,'DefaultAxesFontSize',12)
set(300,'DefaultAxesFontWeight','Bold')
set(300,'DefaultTextfontName','Times')
set(300,'DefaultTextfontSize',12)
set(300,'DefaultTextfontWeight','Bold')
set(300,'Position',[0 0 1000 600])


uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[350 570 300 25],'String','Plot Optical torque');

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'2D-view cutting planes','(x,y)-plane','(x,z)-plane','(y,z)-plane'},...
'Position', [420 500 170 50],...
	  'Callback',{@plottorque,nx,ny,nz,x,y,z,xx,yy,zz,matxytorquex,matxytorquey,matxytorquez,nprint});

if (nprint == 1)
print('-f300','torque2d','-depsc')
end

figure(301)
set(301,'DefaultAxesFontName','Times')
set(301,'DefaultAxesFontSize',12)
set(301,'DefaultAxesFontWeight','Bold')
set(301,'DefaultTextfontName','Times')
set(301,'DefaultTextfontSize',12)
set(301,'DefaultTextfontWeight','Bold')
set(301,'Position',[0 0 1000 600])

quiver3(xx,yy,zz,matxytorquex,matxytorquey,matxytorquez)


xlabel('x')
ylabel('y')
zlabel('z')
title('Density of optical torque')

if (nprint == 1)
print('-f301','torque3d','-depsc')
end

end;  
%%%%%%%%%%%%%%%%%%%%% End  optical torque %%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%% Microscopy %%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nlentille == 1;

if (ntypemic ==0);

if (ntypefile ==0 );
load fourier.mat -ascii
load fourierx.mat -ascii
load fouriery.mat -ascii
load fourierz.mat -ascii
load kxfourier.mat -ascii
elseif (ntypefile ==2 );
kxfourier=h5read(namefileh5,'/Microscopy/kx Fourier');
fourier=h5read(namefileh5,'/Microscopy/Fourier field modulus');
fourierx(:,1)=h5read(namefileh5,'/Microscopy/Fourier field x component real part');
fourierx(:,2)=h5read(namefileh5,'/Microscopy/Fourier field x component imaginary part');
fouriery(:,1)=h5read(namefileh5,'/Microscopy/Fourier field y component real part');
fouriery(:,2)=h5read(namefileh5,'/Microscopy/Fourier field y component imaginary part');
fourierz(:,1)=h5read(namefileh5,'/Microscopy/Fourier field z component real part');
fourierz(:,2)=h5read(namefileh5,'/Microscopy/Fourier field z component imaginary part');
end;

nnfft=max(size(kxfourier));

fourierm=reshape(fourier,nnfft,nnfft);
fourierxc=reshape(fourierx(:,1)+icomp*fourierx(:,2),nnfft,nnfft);
fourieryc=reshape(fouriery(:,1)+icomp*fouriery(:,2),nnfft,nnfft);
fourierzc=reshape(fourierz(:,1)+icomp*fourierz(:,2),nnfft,nnfft);

clear fourier
clear fourierx
clear fouriery
clear fourierz

if (nside==1);

if (ntypefile ==0 );
load fourierinc.mat -ascii
load fourierincx.mat -ascii
load fourierincy.mat -ascii
load fourierincz.mat -ascii
elseif (ntypefile ==2 );
fourierinc=h5read(namefileh5,'/Microscopy/Fourier+incident field modulus');
fourierincx(:,1)=h5read(namefileh5,'/Microscopy/Fourier+incident field x component real part');
fourierincx(:,2)=h5read(namefileh5,'/Microscopy/Fourier+incident field x component imaginary part');
fourierincy(:,1)=h5read(namefileh5,'/Microscopy/Fourier+incident field y component real part');
fourierincy(:,2)=h5read(namefileh5,'/Microscopy/Fourier+incident field y component imaginary part');
fourierincz(:,1)=h5read(namefileh5,'/Microscopy/Fourier+incident field z component real part');
fourierincz(:,2)=h5read(namefileh5,'/Microscopy/Fourier+incident field z component imaginary part');
end;


fourierincm=reshape(fourierinc,nnfft,nnfft);
fourierincxc=reshape(fourierincx(:,1)+icomp*fourierincx(:,2),nnfft,nnfft);
fourierincyc=reshape(fourierincy(:,1)+icomp*fourierincy(:,2),nnfft,nnfft);
fourierinczc=reshape(fourierincz(:,1)+icomp*fourierincz(:,2),nnfft,nnfft);

clear fourierinc
clear fourierincx
clear fourierincy
clear fourierincz

end;

if (ntypefile ==0);
load image.mat -ascii
load imagex.mat -ascii
load imagey.mat -ascii
load imagez.mat -ascii
elseif (ntypefile == 2);
image=h5read(namefileh5,'/Microscopy/Image field modulus');
imagex(:,1)=h5read(namefileh5,'/Microscopy/Image field x component real part');
imagex(:,2)=h5read(namefileh5,'/Microscopy/Image field x component imaginary part');
imagey(:,1)=h5read(namefileh5,'/Microscopy/Image field y component real part');
imagey(:,2)=h5read(namefileh5,'/Microscopy/Image field y component imaginary part');
imagez(:,1)=h5read(namefileh5,'/Microscopy/Image field z component real part');
imagez(:,2)=h5read(namefileh5,'/Microscopy/Image field z component imaginary part');
end;



imagem=reshape(image,nfft,nfft);
imagexc=reshape(imagex(:,1)+icomp*imagex(:,2),nfft,nfft);
imageyc=reshape(imagey(:,1)+icomp*imagey(:,2),nfft,nfft);
imagezc=reshape(imagez(:,1)+icomp*imagez(:,2),nfft,nfft);

clear image
clear imagex
clear imagey
clear imagez

if (nside==1);

if (ntypefile == 0);
load imageinc.mat -ascii
load imageincx.mat -ascii
load imageincy.mat -ascii
load imageincz.mat -ascii
elseif (ntypefile == 2);
imageinc=h5read(namefileh5,'/Microscopy/Image+incident field modulus');
imageincx(:,1)=h5read(namefileh5,'/Microscopy/Image+incident field x component real part');
imageincx(:,2)=h5read(namefileh5,'/Microscopy/Image+incident field x component imaginary part');
imageincy(:,1)=h5read(namefileh5,'/Microscopy/Image+incident field y component real part');
imageincy(:,2)=h5read(namefileh5,'/Microscopy/Image+incident field y component imaginary part');
imageincz(:,1)=h5read(namefileh5,'/Microscopy/Image+incident field z component real part');
imageincz(:,2)=h5read(namefileh5,'/Microscopy/Image+incident field z component imaginary part');
end;

imageincm=reshape(imageinc,nfft,nfft);
imageincxc=reshape(imageincx(:,1)+icomp*imageincx(:,2),nfft,nfft);
imageincyc=reshape(imageincy(:,1)+icomp*imageincy(:,2),nfft,nfft);
imageinczc=reshape(imageincz(:,1)+icomp*imageincz(:,2),nfft,nfft);

clear imageinc
clear imageincx
clear imageincy
clear imageincz
end;

if (ntypefile == 0);
load ximage.mat -ascii
elseif (ntypefile ==2)
ximage=h5read(namefileh5,'/Microscopy/x Image');
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourier %%%%%%%%%%%%%%%%%%%%%%%%%

figure(400)

set(400,'DefaultAxesFontName','Times')
set(400,'DefaultAxesFontSize',12)
set(400,'DefaultAxesFontWeight','Bold')
set(400,'DefaultTextfontName','Times')
set(400,'DefaultTextfontSize',12)
set(400,'DefaultTextfontWeight','Bold')
set(400,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.8 0.8])
  
imagesc(kxfourier,kxfourier,fourierm.^2')
axis xy  
axis equal
axis image

hold on
rectangle('Position',[-numaper -numaper 2*numaper 2*numaper],'Curvature',[1 1],'linewidth',2,'edgecolor','red')

caxis([ 0 max(max(fourierm.^2))])
shading interp
axis equal
axis image
colorbar
  
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)


uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[80 572 400 25],'String','Fourier plane scattered field:')

if (nside == 1);
uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[650 572 80 25],'String','kz>0')
 else  (nside == -1);
uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[650 572 80 25],'String','kz<0')
end;

uicontrol('Style','text','Fontsize',12,'Fontweight','bold','Position',[350 545 200 18],'String','Numerical aperture:')
uicontrol('Style', 'text','Fontsize',12,'Fontweight','bold', 'String', num2str(numaper),'Position', [540 545 40 18]);

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'Intensity','Modulus','x-component','y-component','z-component'},...
'Position', [460 566 150 30],...
	  'Callback',{@plotfourier,numaper,kxfourier,fourierm,fourierxc,fourieryc,fourierzc,nprint});

if (nprint == 1)
print('-f400','fourierpos','-depsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourier +incident %%%%%%%%%%%%%%%%%%%%%%%%%

if (nside==1);

figure(450)

set(450,'DefaultAxesFontName','Times')
set(450,'DefaultAxesFontSize',12)
set(450,'DefaultAxesFontWeight','Bold')
set(450,'DefaultTextfontName','Times')
set(450,'DefaultTextfontSize',12)
set(450,'DefaultTextfontWeight','Bold')
set(450,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.8 0.8])
  
imagesc(kxfourier,kxfourier,fourierincm.^2')
axis xy  
axis equal
axis image

hold on
rectangle('Position',[-numaper -numaper 2*numaper 2*numaper],'Curvature',[1 1],'linewidth',2,'edgecolor','red')

caxis([ 0 max(max(fourierincm^2))])
shading interp
axis equal
axis image
colorbar
  
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)


uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[150 572 300 25],'String','Fourier plane total field:')

uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[650 572 80 25],'String','kz>0')
uicontrol('Style','text','Fontsize',12,'Fontweight','bold','Position',[350 545 200 18],'String','Numerical aperture:')
uicontrol('Style', 'text','Fontsize',12,'Fontweight','bold', 'String', num2str(numaper),'Position', [540 545 40 18]);

uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'Intensity','Modulus','x-component','y-component','z-component'},...
'Position', [460 566 150 30],...
	  'Callback',{@plotfourierinc,numaper,kxfourier,fourierincm,fourierincxc,fourierincyc,fourierinczc,nprint});

if (nprint == 1)
print('-f450','fourierinc','-depsc')
end

end;

%%%%%%%%%%%%%%%%% Image  %%%%%%%%%%%%%%%%%%%%%%


figure(500)

set(500,'DefaultAxesFontName','Times')
set(500,'DefaultAxesFontSize',12)
set(500,'DefaultAxesFontWeight','Bold')
set(500,'DefaultTextfontName','Times')
set(500,'DefaultTextfontSize',12)
set(500,'DefaultTextfontWeight','Bold')
set(500,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.8 0.8])

imagesc(ximage,ximage,imagem.^2')
axis xy
caxis([ 0 max(max(imagem.^2))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)

uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[80 572 400 25],'String','Image plane scattered field:')
if (nside == 1);
uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[650 572 80 25],'String','z>0')
 else  (nside == -1);
uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[650 572 80 25],'String','z<0')
end;
uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'Intensity','Modulus','x-component','y-component','z-component'},...
'Position', [460 561 150 30],...
	  'Callback',{@plotimage,ximage,imagem,imagexc,imageyc,imagezc,nprint});

if (nprint == 1)
print('-f500','image','-depsc')
end


%%%%%%%%%%%%%%%%% Image + incident %%%%%%%%%%%%%%%%%%%%%%

if (nside==1);

figure(550)

set(550,'DefaultAxesFontName','Times')
set(550,'DefaultAxesFontSize',12)
set(550,'DefaultAxesFontWeight','Bold')
set(550,'DefaultTextfontName','Times')
set(550,'DefaultTextfontSize',12)
set(550,'DefaultTextfontWeight','Bold')
set(550,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.8 0.8])

imagesc(ximage,ximage,imageincm.^2')
axis xy
caxis([min(min(imageincm.^2)) max(max(imageincm.^2))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)

uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[150 572 300 25],'String','Image plane total field:')
uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[650 572 80 25],'String','z>0')
uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'Intensity','Modulus','x-component','y-component','z-component'},...
'Position', [460 561 150 30],...
	  'Callback',{@plotimageinc,ximage,imageincm,imageincxc,imageincyc,imageinczc,nprint});

if (nprint == 1)
print('-f550','imageinc','-depsc')
end

end;

elseif (ntypemic == 1);
if (ntypefile == 0);
load imagebf.mat -ascii
load imagebfx.mat -ascii
load imagebfy.mat -ascii
load imagebfz.mat -ascii
elseif(ntypefile == 2);

if (nside==1);
imagebf=h5read(namefileh5,'/Microscopy/Image bright field kz>0 field modulus');
imagebfx(:,1)=h5read(namefileh5,'/Microscopy/Image bright field kz>0 field x component real part');
imagebfx(:,2)=h5read(namefileh5,'/Microscopy/Image bright field kz>0 field x component imaginary part');
imagebfy(:,1)=h5read(namefileh5,'/Microscopy/Image bright field kz>0 field y component real part');
imagebfy(:,2)=h5read(namefileh5,'/Microscopy/Image bright field kz>0 field y component imaginary part');
imagebfz(:,1)=h5read(namefileh5,'/Microscopy/Image bright field kz>0 field z component real part');
imagebfz(:,2)=h5read(namefileh5,'/Microscopy/Image bright field kz>0 field z component imaginary part');
else;
imagebf=h5read(namefileh5,'/Microscopy/Image bright field kz<0 field modulus');
imagebfx(:,1)=h5read(namefileh5,'/Microscopy/Image bright field kz<0 field x component real part');
imagebfx(:,2)=h5read(namefileh5,'/Microscopy/Image bright field kz<0 field x component imaginary part');
imagebfy(:,1)=h5read(namefileh5,'/Microscopy/Image bright field kz<0 field y component real part');
imagebfy(:,2)=h5read(namefileh5,'/Microscopy/Image bright field kz<0 field y component imaginary part');
imagebfz(:,1)=h5read(namefileh5,'/Microscopy/Image bright field kz<0 field z component real part');
imagebfz(:,2)=h5read(namefileh5,'/Microscopy/Image bright field kz<0 field z component imaginary part');
end;

end;


imagem=reshape(imagebf,nfft,nfft);
imagexc=reshape(imagebfx(:,1),nfft,nfft);
imageyc=reshape(imagebfy(:,1),nfft,nfft);
imagezc=reshape(imagebfz(:,1),nfft,nfft);

clear imagebf
clear imagebfx
clear imagebfy
clear imagebfz

if (nside==1);
if (ntypefile == 0);
load imageincbf.mat -ascii
load imageincbfx.mat -ascii
load imageincbfy.mat -ascii
load imageincbfz.mat -ascii
elseif (ntypefile == 2);
imageincbf=h5read(namefileh5,'/Microscopy/Image+incident bright field kz>0 field modulus');
imageincbfx(:,1)=h5read(namefileh5,'/Microscopy/Image+incident bright field kz>0 field x component real part');
imageincbfx(:,2)=h5read(namefileh5,'/Microscopy/Image+incident bright field kz>0 field x component imaginary part');
imageincbfy(:,1)=h5read(namefileh5,'/Microscopy/Image+incident bright field kz>0 field y component real part');
imageincbfy(:,2)=h5read(namefileh5,'/Microscopy/Image+incident bright field kz>0 field y component imaginary part');
imageincbfz(:,1)=h5read(namefileh5,'/Microscopy/Image+incident bright field kz>0 field z component real part');
imageincbfz(:,2)=h5read(namefileh5,'/Microscopy/Image+incident bright field kz>0 field z component imaginary part');
end;

imageincm=reshape(imageincbf,nfft,nfft);
imageincxc=reshape(imageincbfx(:,1),nfft,nfft);
imageincyc=reshape(imageincbfy(:,1),nfft,nfft);
imageinczc=reshape(imageincbfz(:,1),nfft,nfft);

clear imageincbf
clear imageincbfx
clear imageincbfy
clear imageincbfz

end;

if (ntypefile == 0);
load ximage.mat -ascii
elseif (ntypefile ==2)
ximage=h5read(namefileh5,'/Microscopy/x Image');
end;

%%%%%%%%%%%%%%%%% Image  %%%%%%%%%%%%%%%%%%%%%%


figure(500)

set(500,'DefaultAxesFontName','Times')
set(500,'DefaultAxesFontSize',12)
set(500,'DefaultAxesFontWeight','Bold')
set(500,'DefaultTextfontName','Times')
set(500,'DefaultTextfontSize',12)
set(500,'DefaultTextfontWeight','Bold')
set(500,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.75 0.75])

imagesc(ximage,ximage,imagem.^2')
axis xy
caxis([ 0 max(max(imagem.^2))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)

uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[80 572 400 25],'String','Image plane scattered field:')
if (nside == 1);
uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[650 572 80 25],'String','z>0')
 else  (nside == -1);
uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[650 572 80 25],'String','z<0')
end;
uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'Intensity','Modulus','x-component','y-component','z-component'},...
'Position', [460 561 150 30],...
	  'Callback',{@plotimagereal,ximage,imagem,imagexc,imageyc,imagezc,nprint});

if (nprint == 1)
print('-f500','imagewf','-depsc')
end

%%%%%%%%%%%%%%%%% Image + incident %%%%%%%%%%%%%%%%%%%%%%

if (nside==1);

figure(550)

set(550,'DefaultAxesFontName','Times')
set(550,'DefaultAxesFontSize',12)
set(550,'DefaultAxesFontWeight','Bold')
set(550,'DefaultTextfontName','Times')
set(550,'DefaultTextfontSize',12)
set(550,'DefaultTextfontWeight','Bold')
set(550,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.75 0.75])

imagesc(ximage,ximage,imageincm.^2')
axis xy
caxis([min(min(imageincm.^2)) max(max(imageincm.^2))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)

uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[150 572 300 25],'String','Image plane total field:')
uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[650 572 80 25],'String','z>0')
uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'Intensity','Modulus','x-component','y-component','z-component'},...
'Position', [460 561 150 30],...
	  'Callback',{@plotimageincreal,ximage,imageincm,imageincxc,imageincyc,imageinczc,nprint});

if (nprint == 1)
print('-f550','imageincwf','-depsc')
end

end;


figure(560)
set(560,'DefaultAxesFontName','Times')
set(560,'DefaultAxesFontSize',12)
set(560,'DefaultAxesFontWeight','Bold')
set(560,'DefaultTextfontName','Times')
set(560,'DefaultTextfontSize',12)
set(560,'DefaultTextfontWeight','Bold')
set(560,'Position',[0 0 1000 600])

if (ntypefile == 0);
load kxincidentbf.mat -ascii
load kyincidentbf.mat -ascii
elseif (ntypefile ==2)
kxincidentbf=h5read(namefileh5,'/Microscopy/kx incident bf');
kyincidentbf=h5read(namefileh5,'/Microscopy/ky incident bf');
end;

plot(kxincidentbf,kyincidentbf,'b+')
title('Position in Fourier space of all the incident field to compute bf')
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)
axis xy  
axis equal
axis image

hold on
rectangle('Position',[-numaperinc -numaperinc 2*numaperinc 2*numaperinc],'Curvature',[1 1],'linewidth',2,'edgecolor','red')

if (nprint == 1)
print('-f560','angleincbf','-depsc')
end


elseif (ntypemic == 2);

if (ntypefile == 0);
load imagedf.mat -ascii
load imagedfx.mat -ascii
load imagedfy.mat -ascii
load imagedfz.mat -ascii
elseif(ntypefile == 2);

if (nside==1);
imagedf=h5read(namefileh5,'/Microscopy/Image dark field kz>0 field modulus');
imagedfx(:,1)=h5read(namefileh5,'/Microscopy/Image dark field kz>0 field x component real part');
imagedfx(:,2)=h5read(namefileh5,'/Microscopy/Image dark field kz>0 field x component imaginary part');
imagedfy(:,1)=h5read(namefileh5,'/Microscopy/Image dark field kz>0 field y component real part');
imagedfy(:,2)=h5read(namefileh5,'/Microscopy/Image dark field kz>0 field y component imaginary part');
imagedfz(:,1)=h5read(namefileh5,'/Microscopy/Image dark field kz>0 field z component real part');
imagedfz(:,2)=h5read(namefileh5,'/Microscopy/Image dark field kz>0 field z component imaginary part');
else;
imagedf=h5read(namefileh5,'/Microscopy/Image dark field kz<0 field modulus');
imagedfx(:,1)=h5read(namefileh5,'/Microscopy/Image dark field kz<0 field x component real part');
imagedfx(:,2)=h5read(namefileh5,'/Microscopy/Image dark field kz<0 field x component imaginary part');
imagedfy(:,1)=h5read(namefileh5,'/Microscopy/Image dark field kz<0 field y component real part');
imagedfy(:,2)=h5read(namefileh5,'/Microscopy/Image dark field kz<0 field y component imaginary part');
imagedfz(:,1)=h5read(namefileh5,'/Microscopy/Image dark field kz<0 field z component real part');
imagedfz(:,2)=h5read(namefileh5,'/Microscopy/Image dark field kz<0 field z component imaginary part');
end;

end;





imagem=reshape(imagedf,nfft,nfft);
imagexc=reshape(imagedfx(:,1),nfft,nfft);
imageyc=reshape(imagedfy(:,1),nfft,nfft);
imagezc=reshape(imagedfz(:,1),nfft,nfft);

clear imagedf
clear imagedfx
clear imagedfy
clear imagedfz

if (nside==1);
if (ntypefile == 0);
load imageincdf.mat -ascii
load imageincdfx.mat -ascii
load imageincdfy.mat -ascii
load imageincdfz.mat -ascii
elseif (ntypefile == 2);
imageincdf=h5read(namefileh5,'/Microscopy/Image+incident dark field kz>0 field modulus');
imageincdfx(:,1)=h5read(namefileh5,'/Microscopy/Image+incident dark field kz>0 field x component real part');
imageincdfx(:,2)=h5read(namefileh5,'/Microscopy/Image+incident dark field kz>0 field x component imaginary part');
imageincdfy(:,1)=h5read(namefileh5,'/Microscopy/Image+incident dark field kz>0 field y component real part');
imageincdfy(:,2)=h5read(namefileh5,'/Microscopy/Image+incident dark field kz>0 field y component imaginary part');
imageincdfz(:,1)=h5read(namefileh5,'/Microscopy/Image+incident dark field kz>0 field z component real part');
imageincdfz(:,2)=h5read(namefileh5,'/Microscopy/Image+incident dark field kz>0 field z component imaginary part');
end;

imageincm=reshape(imageincdf,nfft,nfft);
imageincxc=reshape(imageincdfx(:,1),nfft,nfft);
imageincyc=reshape(imageincdfy(:,1),nfft,nfft);
imageinczc=reshape(imageincdfz(:,1),nfft,nfft);

end;

clear imageincdf
clear imageincdfx
clear imageincdfy
clear imageincdfz

if (ntypefile == 0);
load ximage.mat -ascii
elseif (ntypefile ==2)
ximage=h5read(namefileh5,'/Microscopy/x Image');
end;

%%%%%%%%%%%%%%%%% Image  %%%%%%%%%%%%%%%%%%%%%%


figure(500)

set(500,'DefaultAxesFontName','Times')
set(500,'DefaultAxesFontSize',12)
set(500,'DefaultAxesFontWeight','Bold')
set(500,'DefaultTextfontName','Times')
set(500,'DefaultTextfontSize',12)
set(500,'DefaultTextfontWeight','Bold')
set(500,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.75 0.75])

imagesc(ximage,ximage,imagem.^2')
axis xy
caxis([ 0 max(max(imagem.^2))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)

uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[80 572 400 25],'String','Image plane scattered field:')
if (nside == 1);
uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[650 572 80 25],'String','z>0')
 else  (nside == -1);
uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[650 572 80 25],'String','z<0')
end;
uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'Intensity','Modulus','x-component','y-component','z-component'},...
'Position', [460 561 150 30],...
	  'Callback',{@plotimagereal,ximage,imagem,imagexc,imageyc,imagezc,nprint});

if (nprint == 1)
print('-f500','imagewf','-depsc')
end


%%%%%%%%%%%%%%%%% Image + incident %%%%%%%%%%%%%%%%%%%%%%

if (nside==1);

figure(550)

set(550,'DefaultAxesFontName','Times')
set(550,'DefaultAxesFontSize',12)
set(550,'DefaultAxesFontWeight','Bold')
set(550,'DefaultTextfontName','Times')
set(550,'DefaultTextfontSize',12)
set(550,'DefaultTextfontWeight','Bold')
set(550,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.75 0.75])

imagesc(ximage,ximage,imageincm.^2')
axis xy
caxis([min(min(imageincm.^2)) max(max(imageincm.^2))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)

uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[150 572 300 25],'String','Image plane total field:')
uicontrol('Style','text','Fontsize',16,'Fontweight','bold','Position',[650 572 80 25],'String','z>0')
uicontrol('Style', 'popupmenu','Fontsize',12,'String',...
{'Intensity','Modulus','x-component','y-component','z-component'},...
'Position', [460 561 150 30],...
	  'Callback',{@plotimageincreal,ximage,imageincm,imageincxc,imageincyc,imageinczc,nprint});
	
if (nprint == 1)
print('-f550','imageincwf','-depsc')
end

end;



figure(560)
set(560,'DefaultAxesFontName','Times')
set(560,'DefaultAxesFontSize',12)
set(560,'DefaultAxesFontWeight','Bold')
set(560,'DefaultTextfontName','Times')
set(560,'DefaultTextfontSize',12)
set(560,'DefaultTextfontWeight','Bold')
set(560,'Position',[0 0 1000 600])

if (ntypefile == 0);
load kxincidentdf.mat -ascii
load kyincidentdf.mat -ascii
elseif (ntypefile ==2)
kxincidentdf=h5read(namefileh5,'/Microscopy/kx incident df');
kyincidentdf=h5read(namefileh5,'/Microscopy/ky incident df');
end;

plot(kxincidentdf,kyincidentdf,'o','MarkerSize',5,...
     'MarkerEdgeColor','blue','MarkerFaceColor',[0 0 1])
title('Position in Fourier space of all the incident field to compute df')
xlabel('$k_x/k_0$','Interpreter','latex','Fontsize',18)
ylabel('$k_y/k_0$','Interpreter','latex','Fontsize',18)
axis xy  
axis equal
axis image

hold on
rectangle('Position',[-numaperinc -numaperinc 2*numaperinc 2*numaperinc],'Curvature',[1 1],'linewidth',1,'edgecolor','red')

if (nprint == 1)
print('-f560','angleincdf','-depsc')
end

end;
	
end;
