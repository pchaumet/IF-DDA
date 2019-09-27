      subroutine incidentarbitrary(xs,ys,zs,aretecube,FF0,nxm,nym,nzm
     $     ,nbsphere,nstop,namefile,infostr)
      implicit none
      integer ii,iii,i,nx,ny,nz,nxy,nxm,nym,nzm,nbsphere,imin,jmin
     $     ,kmin,nstop,ierror
      double precision dx,dy,dz,xmin,ymin,zmin,xmax,ymax,zmax,exr,exi
     $     ,eyr,eyi,ezr,ezi,xd,yd,zd,aretecube,dmax

      double complex Ex(8),Ey(8),Ez(8)

      double precision xs(nxm*nym*nzm),ys(nxm*nym*nzm),zs(nxm*nym*nzm)
      double complex FF0(3*nxm*nym*nzm),uncomp,icomp
      character(64) infostr
      character(64) namefile

      uncomp=(1.d0,0.d0)
      icomp=(0.d0,1.d0)

c     read the input file
      open(15,file=namefile,status='old',iostat=ierror) 
      if (ierror.ne.0) then
         infostr='arbitrary incident: name file does not exist'
         nstop=1
         return
      endif

      read(15,*) nx,ny,nz           
      read(15,*) dx,dy,dz
      read(15,*) xmin,ymin,zmin
      close(15)
      nxy=nx*ny

      dx=dx*1.d-9
      dy=dy*1.d-9
      dz=dz*1.d-9
      dmax=max(dx,dy,dz)
      if (aretecube.le.dmax) then
         infostr='object mesh size < mesh size incident field'
         nstop=1
         return
      endif


      xmin=xmin*1.d-9
      ymin=ymin*1.d-9
      zmin=zmin*1.d-9

      xmax=xmin+dble(nx-1)*dx
      ymax=ymin+dble(ny-1)*dy
      zmax=zmin+dble(nz-1)*dz

c     verifier que l'objet est dans la boite

      if (xs(1).lt.xmin.or.xs(nbsphere).gt.xmax .or.
     $     ys(1).lt.ymin.or.ys(nbsphere).gt.ymax .or.
     $     zs(1).lt.zmin.or.zs(nbsphere).gt.zmax) then
         infostr='objet not inside the incident box'
         nstop=1
         write(*,*) 'objet not inside the incident box'
         write(*,*)  'position min and max',xmin,xmax,ymin ,ymax,zmin
     $        ,zmax
         write(*,*) 'position of the corners',xs(1),xs(nbsphere),ys(1)
     $        ,ys(nbsphere),zs(1),zs(nbsphere)
         return
      endif
c     ecriture des fichiers sous cette structure
c     do i=1,nz
c     do j=1,ny
c     do k=1,nx
c     x=x0+dble(k-1)*dx
c     y=y0+dble(j-1)*dy
c     z=z0+dble(i-1)*dz

      open(11, file='Exr.mat', status='old', form='formatted', access
     $     ='direct', recl=22,iostat=ierror) 
      if (ierror.ne.0) then
         infostr='arbitrary incident: Exr.mat does not exist'
         nstop=1
         return
      endif
      open(12, file='Exi.mat', status='old', form='formatted', access
     $     ='direct', recl=22,iostat=ierror)
      if (ierror.ne.0) then
         infostr='arbitrary incident: Exi.mat does not exist'
         nstop=1
         return
      endif
      open(13, file='Eyr.mat', status='old', form='formatted', access
     $     ='direct', recl=22,iostat=ierror)
      if (ierror.ne.0) then
         infostr='arbitrary incident: Eyr.mat does not exist'
         nstop=1
         return
      endif
      open(14, file='Eyi.mat', status='old', form='formatted', access
     $     ='direct', recl=22,iostat=ierror)
      if (ierror.ne.0) then
         infostr='arbitrary incident: Eyi.mat does not exist'
         nstop=1
         return
      endif
      open(15, file='Ezr.mat', status='old', form='formatted', access
     $     ='direct', recl=22,iostat=ierror)
      if (ierror.ne.0) then
         infostr='arbitrary incident: Ezr.mat does not exist'
         nstop=1
         return
      endif
      open(16, file='Ezi.mat', status='old', form='formatted', access
     $     ='direct', recl=22,iostat=ierror)
      if (ierror.ne.0) then
         infostr='arbitrary incident: Ezi.mat does not exist'
         nstop=1
         return
      endif
      do i=1,nbsphere
         iii=3*(i-1)
         imin=IDNINT(dint((xs(i)-xmin)/dx))+1
         jmin=IDNINT(dint((ys(i)-ymin)/dy))+1
         kmin=IDNINT(dint((zs(i)-zmin)/dz))+1

         xd=(xs(i)-xmin)/dx-dble(imin-1)
         yd=(ys(i)-ymin)/dy-dble(jmin-1)
         zd=(zs(i)-zmin)/dz-dble(kmin-1)
         
         if (xd.le.dx/1000.d0.and.yd.le.dy/1000.d0.and.zd.le.dz/1000.d0)
     $        then
c     on tombe sur une maille
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(1)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(1)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(1)=Ezr*uncomp+icomp*Ezi            
            FF0(iii+1)=Ex(1)
            FF0(iii+2)=Ey(1)
            FF0(iii+3)=Ez(1)
         elseif (xd.le.dx/1000.d0.and.yd.le.dy/1000.d0) then

            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(1)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(1)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(1)=Ezr*uncomp+icomp*Ezi   

            ii=imin+nx*(jmin-1)+nxy*kmin
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(5)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(5)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(5)=Ezr*uncomp+icomp*Ezi

            FF0(iii+1)=(1.d0-zd)*Ex(1)+zd*Ex(5)
            FF0(iii+2)=(1.d0-zd)*Ey(1)+zd*Ey(5)
            FF0(iii+3)=(1.d0-zd)*Ez(1)+zd*Ez(5)

         elseif (xd.le.dx/1000.d0.and.zd.le.dz/1000.d0) then
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(1)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(1)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(1)=Ezr*uncomp+icomp*Ezi   

            ii=imin+nx*jmin+nxy*(kmin-1)
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(3)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(3)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(3)=Ezr*uncomp+icomp*Ezi

            FF0(iii+1)=(1.d0-yd)*Ex(1)+yd*Ex(3)
            FF0(iii+2)=(1.d0-yd)*Ey(1)+yd*Ey(3)
            FF0(iii+3)=(1.d0-yd)*Ez(1)+yd*Ez(3)


         elseif (yd.le.dy/1000.d0.and.zd.le.dz/1000.d0) then
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(1)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(1)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(1)=Ezr*uncomp+icomp*Ezi   

            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(2)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(2)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(2)=Ezr*uncomp+icomp*Ezi

            FF0(iii+1)=Ex(1)*(1.d0-xd)+Ex(2)*xd
            FF0(iii+2)=Ey(1)*(1.d0-xd)+Ey(2)*xd
            FF0(iii+3)=Ez(1)*(1.d0-xd)+Ez(2)*xd

         elseif (xd.le.dx/1000.d0) then


c     1 coin
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(1)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(1)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(1)=Ezr*uncomp+icomp*Ezi

c     3 coin
            ii=imin+nx*jmin+nxy*(kmin-1)
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(3)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(3)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(3)=Ezr*uncomp+icomp*Ezi

c     5 coin
            ii=imin+nx*(jmin-1)+nxy*kmin
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(5)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(5)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(5)=Ezr*uncomp+icomp*Ezi


c     7 coin
            ii=imin+nx*jmin+nxy*kmin
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(7)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(7)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(7)=Ezr*uncomp+icomp*Ezi


            FF0(iii+1)=(1.d0-zd)*((1-yd)*Ex(1)+yd*Ex(3))+zd*(Ex(5)
     $           *(1.d0-yd)+Ex(7)*yd)
            FF0(iii+2)=(1.d0-zd)*((1-yd)*Ey(1)+yd*Ey(3))+zd*(Ey(5)
     $           *(1.d0-yd)+Ey(7)*yd)
            FF0(iii+3)=(1.d0-zd)*((1-yd)*Ez(1)+yd*Ez(3))+zd*(Ez(5)
     $           *(1.d0-yd)+Ez(7)*yd)

         elseif (yd.le.dy/1000.d0) then

            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(1)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(1)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(1)=Ezr*uncomp+icomp*Ezi

            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(2)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(2)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(2)=Ezr*uncomp+icomp*Ezi

            ii=imin+nx*(jmin-1)+nxy*kmin
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(5)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(5)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(5)=Ezr*uncomp+icomp*Ezi

            ii=imin+1+nx*(jmin-1)+nxy*kmin
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(6)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(6)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(6)=Ezr*uncomp+icomp*Ezi

            FF0(iii+1)=(1.d0-zd)*((Ex(1)*(1.d0-xd)+Ex(2)*xd)) +zd
     $           *((Ex(5)*(1.d0-xd)+Ex(6)*xd))
            FF0(iii+2)=(1.d0-zd)*((Ey(1)*(1.d0-xd)+Ey(2)*xd)) +zd
     $           *((Ey(5)*(1.d0-xd)+Ey(6)*xd))
            FF0(iii+3)=(1.d0-zd)*((Ez(1)*(1.d0-xd)+Ez(2)*xd)) +zd
     $           *((Ez(5)*(1.d0-xd)+Ez(6)*xd))

         elseif (zd.le.dz/1000.d0) then
c     1 coin
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(1)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(1)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(1)=Ezr*uncomp+icomp*Ezi

c     2 coin
            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(2)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(2)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(2)=Ezr*uncomp+icomp*Ezi

c     3 coin
            ii=imin+nx*jmin+nxy*(kmin-1)
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(3)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(3)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(3)=Ezr*uncomp+icomp*Ezi

c     4 coin
            ii=imin+1+nx*jmin+nxy*(kmin-1)
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(4)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(4)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(4)=Ezr*uncomp+icomp*Ezi

            FF0(iii+1)=(1.d0-yd)*(Ex(1)*(1.d0-xd)+Ex(2)*xd)+yd*(Ex(3)
     $           *(1.d0-xd) + Ex(4)*xd)
            FF0(iii+2)=(1.d0-yd)*(Ey(1)*(1.d0-xd)+Ey(2)*xd)+yd*(Ey(3)
     $           *(1.d0-xd)+Ey(4)*xd)
            FF0(iii+3)=(1.d0-yd)*(Ez(1)*(1.d0-xd)+Ez(2)*xd)+yd*(Ez(3)
     $           *(1.d0-xd)+Ez(4)*xd)

         else


c     1 coin
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)

c     WRITE(*,*) 'rr',(xs(i)-xmin)/dx,imin,jmin,kmin,ii

            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(1)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(1)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(1)=Ezr*uncomp+icomp*Ezi

c     2 coin
            
            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
c     WRITE(*,*) 'rr2',ii
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(2)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(2)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(2)=Ezr*uncomp+icomp*Ezi

c     3 coin
            ii=imin+nx*jmin+nxy*(kmin-1)
c     WRITE(*,*) 'rr3',ii
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(3)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(3)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(3)=Ezr*uncomp+icomp*Ezi

c     4 coin
            ii=imin+1+nx*jmin+nxy*(kmin-1)
c     WRITE(*,*) 'rr4',ii
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(4)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(4)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(4)=Ezr*uncomp+icomp*Ezi

c     5 coin
            ii=imin+nx*(jmin-1)+nxy*kmin
c     WRITE(*,*) 'rr5',ii
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(5)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(5)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(5)=Ezr*uncomp+icomp*Ezi

c     6 coin
            ii=imin+1+nx*(jmin-1)+nxy*kmin
c     WRITE(*,*) 'rr6',ii
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(6)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(6)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(6)=Ezr*uncomp+icomp*Ezi

c     7 coin
            ii=imin+nx*jmin+nxy*kmin
c     WRITE(*,*) 'rr7',ii
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(7)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(7)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(7)=Ezr*uncomp+icomp*Ezi

c     8 coin
            ii=imin+1+nx*jmin+nxy*kmin
c     WRITE(*,*) 'rr8',ii
            read(11,FMT='(D22.15)',rec=ii) Exr
            read(12,FMT='(D22.15)',rec=ii) Exi
            Ex(8)=Exr*uncomp+icomp*Exi
            read(13,FMT='(D22.15)',rec=ii) Eyr
            read(14,FMT='(D22.15)',rec=ii) Eyi
            Ey(8)=Eyr*uncomp+icomp*Eyi
            read(15,FMT='(D22.15)',rec=ii) Ezr
            read(16,FMT='(D22.15)',rec=ii) Ezi
            Ez(8)=Ezr*uncomp+icomp*Ezi

c     interpolation trilineaire      
c     c_{00} = V[x_0,y_0, z_0] (1 - x_d) + V[x_1, y_0, z_0] x_d 
c     c_{10} = V[x_0,y_1, z_0] (1 - x_d) + V[x_1, y_1, z_0] x_d 
c     c_{01} = V[x_0,y_0, z_1] (1 - x_d) + V[x_1, y_0, z_1] x_d 
c     c_{11} = V[x_0,y_1, z_1] (1 - x_d) + V[x_1, y_1, z_1] x_d 
c     c_0 = c_{00}(1 - y_d) + c_{10}y_d
c     c_1 = c_{01}(1 - y_d) + c_{11}y_d
c     c = c_0(1 - z_d) + c_1z_d .    

            FF0(iii+1)=(1.d0-zd)*((1-yd)*(Ex(1)*(1.d0-xd) +Ex(2)*xd) +yd
     $           *(Ex(3)*(1.d0-xd) + Ex(4)*xd )) +zd*((Ex(5)*(1.d0-xd) +
     $           Ex(6) *xd)*(1.d0-yd) + (Ex(7)*(1.d0-xd) + Ex(8)*xd)*yd)

            FF0(iii+2)=(1.d0-zd)*((1-yd)*(Ey(1)*(1.d0-xd) +Ey(2)*xd) +yd
     $           *(Ey(3)*(1.d0-xd) + Ey(4)*xd )) +zd*((Ey(5)*(1.d0-xd) +
     $           Ey(6) *xd)*(1.d0-yd) + (Ey(7)*(1.d0-xd) + Ey(8)*xd)*yd)

            FF0(iii+3)=(1.d0-zd)*((1-yd)*(Ez(1)*(1.d0-xd) +Ez(2)*xd) +yd
     $           *(Ez(3)*(1.d0-xd) + Ez(4)*xd )) +zd*((Ez(5)*(1.d0-xd) +
     $           Ez(6) *xd)*(1.d0-yd) + (Ez(7)*(1.d0-xd) + Ez(8)*xd)*yd)

         endif

      enddo
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      end
c     *************************************************************
c     *************************************************************
c     *************************************************************
c     subroutine qui ecrit calculs les derivees d'apres le fichier
c     d'entree du champ
      subroutine incidentarbitrarydercreate(namefile)
      implicit none
      integer i,j,k,ii,nx,ny,nz,nxy,ierror
      double precision dx,dy,dz,Exr,Exi,xmin,ymin,zmin
      double complex icomp,uncomp,FF(1000),FFder(1000)
      character(64) namefile

      uncomp=(1.d0,0.d0)
      icomp=(0.d0,1.d0)

c     read the input file
      open(15,file=namefile,status='old')     
      read(15,*) nx,ny,nz           
      read(15,*) dx,dy,dz
      read(15,*) xmin,ymin,zmin
      close(15)
      nxy=nx*ny

      dx=dx*1.d-9
      dy=dy*1.d-9
      dz=dz*1.d-9
      xmin=xmin*1.d-9
      ymin=ymin*1.d-9
      zmin=zmin*1.d-9

c      write(*,*) 'derive',nx,ny,nz,dx,dy,dz,xmin,ymin,zmin
      open(11, file='Exr.mat', status='old', form='formatted', access
     $     ='direct', recl=22)
      open(12, file='Exi.mat', status='old', form='formatted', access
     $     ='direct', recl=22)
      open(13, file='Eyr.mat', status='old', form='formatted', access
     $     ='direct', recl=22)
      open(14, file='Eyi.mat', status='old', form='formatted', access
     $     ='direct', recl=22)
      open(15, file='Ezr.mat', status='old', form='formatted', access
     $     ='direct', recl=22)
      open(16, file='Ezi.mat', status='old', form='formatted', access
     $     ='direct', recl=22)

      open(21, file='Exdx.mat', status='new', form='unformatted', access
     $     ='direct', recl=16,iostat=ierror)   
C     TEST SI LE FICHIER EXISTE DEJA
      if (ierror.ne.0) return

      open(22, file='Eydx.mat', status='new', form='unformatted', access
     $     ='direct', recl=16)
      open(23, file='Ezdx.mat', status='new', form='unformatted', access
     $     ='direct', recl=16)
      open(24, file='Exdy.mat', status='new', form='unformatted', access
     $     ='direct', recl=16)   
      open(25, file='Eydy.mat', status='new', form='unformatted', access
     $     ='direct', recl=16)
      open(26, file='Ezdy.mat', status='new', form='unformatted', access
     $     ='direct', recl=16)
      open(27, file='Exdz.mat', status='new', form='unformatted', access
     $     ='direct', recl=16)   
      open(28, file='Eydz.mat', status='new', form='unformatted', access
     $     ='direct', recl=16)
      open(29, file='Ezdz.mat', status='new', form='unformatted', access
     $     ='direct', recl=16)

c     calcul derivee suivant x
      do i=1,nz
         do j=1,ny
            do k=1,nx
               ii=k+nx*(j-1)+nxy*(i-1)
               read(11,FMT='(D22.15)',rec=ii) Exr
               read(12,FMT='(D22.15)',rec=ii) Exi
               FF(k)=Exr+icomp*Exi
            enddo

            call deriveFF(FF,FFder,nx,dx)
c     ecriture de Exdx
            do k=1,nx
               ii=k+nx*(j-1)+nxy*(i-1)
               write(21,rec=ii) FFder(k)            
            enddo
            do k=1,nx
               ii=k+nx*(j-1)+nxy*(i-1)
               read(13,FMT='(D22.15)',rec=ii) Exr
               read(14,FMT='(D22.15)',rec=ii) Exi
               FF(k)=Exr+icomp*Exi
            enddo
            call deriveFF(FF,FFder,nx,dx)
c     ecriture de Eydx
             do k=1,nx
               ii=k+nx*(j-1)+nxy*(i-1)
               write(22,rec=ii) FFder(k)
            enddo

            do k=1,nx
               ii=k+nx*(j-1)+nxy*(i-1)
               read(15,FMT='(D22.15)',rec=ii) Exr
               read(16,FMT='(D22.15)',rec=ii) Exi
               FF(k)=Exr+icomp*Exi
            enddo
            call deriveFF(FF,FFder,nx,dx)
c     ecriture de Exdx
             do k=1,nx
               ii=k+nx*(j-1)+nxy*(i-1)
               write(23,rec=ii) FFder(k)
            enddo

         enddo
      enddo

c     calcul des derivees suivant y
      do i=1,nz
         do k=1,nx
            do j=1,ny          
               ii=k+nx*(j-1)+nxy*(i-1)
               read(11,FMT='(D22.15)',rec=ii) Exr
               read(12,FMT='(D22.15)',rec=ii) Exi
               FF(j)=Exr+icomp*Exi
            enddo
            call deriveFF(FF,FFder,ny,dy)
c     ecriture de Exdy
             do j=1,ny
               ii=k+nx*(j-1)+nxy*(i-1)
               write(24,rec=ii) FFder(j)
            enddo

            do j=1,ny          
               ii=k+nx*(j-1)+nxy*(i-1)
               read(13,FMT='(D22.15)',rec=ii) Exr
               read(14,FMT='(D22.15)',rec=ii) Exi
               FF(j)=Exr+icomp*Exi
            enddo
            call deriveFF(FF,FFder,ny,dy)
c     ecriture de Eydy
             do j=1,ny
               ii=k+nx*(j-1)+nxy*(i-1)
               write(25,rec=ii) FFder(j)
            enddo
            do j=1,ny          
               ii=k+nx*(j-1)+nxy*(i-1)
               read(15,FMT='(D22.15)',rec=ii) Exr
               read(16,FMT='(D22.15)',rec=ii) Exi
               FF(j)=Exr+icomp*Exi
            enddo
            call deriveFF(FF,FFder,ny,dy)
c     ecriture de Ezdy
             do j=1,ny
               ii=k+nx*(j-1)+nxy*(i-1)
               write(26,rec=ii) FFder(j)
            enddo
         enddo
      enddo

c     calcul des derivees suivant z
    
      do k=1,nx
         do j=1,ny 
            do i=1,nz         
               ii=k+nx*(j-1)+nxy*(i-1)
               read(11,FMT='(D22.15)',rec=ii) Exr
               read(12,FMT='(D22.15)',rec=ii) Exi
               FF(i)=Exr+icomp*Exi
            enddo
            call deriveFF(FF,FFder,nz,dz)
c     ecriture de Exdz
             do i=1,nz
               ii=k+nx*(j-1)+nxy*(i-1)
               write(27,rec=ii) FFder(i)
c               write(*,*) 'derz',i,FF(i),FFder(i),k,j
            enddo

            do i=1,nz         
               ii=k+nx*(j-1)+nxy*(i-1)
               read(13,FMT='(D22.15)',rec=ii) Exr
               read(14,FMT='(D22.15)',rec=ii) Exi
               FF(i)=Exr+icomp*Exi
            enddo
            call deriveFF(FF,FFder,nz,dz)
c     ecriture de Eydz
             do i=1,nz
               ii=k+nx*(j-1)+nxy*(i-1)
               write(28,rec=ii) FFder(i)
            enddo
            do i=1,nz          
               ii=k+nx*(j-1)+nxy*(i-1)
               read(15,FMT='(D22.15)',rec=ii) Exr
               read(16,FMT='(D22.15)',rec=ii) Exi
               FF(i)=Exr+icomp*Exi
            enddo
            call deriveFF(FF,FFder,nz,dz)
c     ecriture de Ezdz
            do i=1,nz
               ii=k+nx*(j-1)+nxy*(i-1)
               write(29,rec=ii) FFder(i)
            enddo
         enddo
      enddo

      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)
      close(29)

      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      end
c     ****************************************************************
c     ****************************************************************
c     ****************************************************************
      subroutine deriveFF(FF,FFder,nx,dx)
      implicit none
      integer i,nx,nordre,nderiv
      double precision dx,pos(0:4),xderiv
      double complex FF(1000),FFder(1000),f(0:4),deriv

      nordre=5
      nderiv=1
      pos(0)=0.d0
      pos(1)=dx
      pos(2)=dx*2.d0
      pos(3)=dx*3.d0
      pos(4)=dx*4.d0
      
      f(0)=FF(1)
      f(1)=FF(2)
      f(2)=FF(3)
      f(3)=FF(4)
      f(4)=FF(5)
      xderiv=pos(0)
      call derivative(nordre,nderiv,f,pos,dx,xderiv,deriv)
      FFder(1)=deriv

      xderiv=pos(1)
      call derivative(nordre,nderiv,f,pos,dx,xderiv,deriv)
      FFder(2)=deriv

      f(0)=FF(nx-4)
      f(1)=FF(nx-3)
      f(2)=FF(nx-2)
      f(3)=FF(nx-1)
      f(4)=FF(nx)
      xderiv=pos(4)
      call derivative(nordre,nderiv,f,pos,dx,xderiv,deriv)
      FFder(nx)=deriv

      xderiv=pos(3)
      call derivative(nordre,nderiv,f,pos,dx,xderiv,deriv)
      FFder(nx-1)=deriv

      do i=3,nx-2
         f(0)=FF(i-2)
         f(1)=FF(i-1)
         f(2)=FF(i)
         f(3)=FF(i+1)
         f(4)=FF(i+2)
         xderiv=pos(2)
         call derivative(nordre,nderiv,f,pos,dx,xderiv,deriv)
         FFder(i)=deriv         
      enddo

      end
c     *********************************************************
c     *********************************************************
c     *********************************************************
      subroutine incidentarbitraryder1(test,x,y,z,namefile,Ederint)
      implicit none
      integer test,nx,ny,nz,nxy,ii,jj,kk,imin,jmin,kmin
      double precision x,y,z,dx,dy,dz,xmin,ymin,zmin,xmax,ymax,zmax,xd
     $     ,yd,zd
      double complex Ederint(3,3),Eder(8,3,3)    
      character(64) namefile

c     read the input file
      open(15,file=namefile,status='old') 
      read(15,*) nx,ny,nz           
      read(15,*) dx,dy,dz
      read(15,*) xmin,ymin,zmin
      close(15)
      nxy=nx*ny

      dx=dx*1.d-9
      dy=dy*1.d-9
      dz=dz*1.d-9
      xmin=xmin*1.d-9
      ymin=ymin*1.d-9
      zmin=zmin*1.d-9

      xmax=xmin+dble(nx-1)*dx
      ymax=ymin+dble(ny-1)*dy
      zmax=zmin+dble(nz-1)*dz

c     verifier que l'objet est dans la boite

      if (x.lt.xmin.or.x.gt.xmax .or. y.lt.ymin.or.y.gt.ymax .or.
     $     z.lt.zmin.or.z.gt.zmax) then
         write(*,*) 'objet not inside the incident box'
         write(*,*)  'position min and max',xmin,xmax,ymin ,ymax,zmin
     $        ,zmax
         write(*,*) 'coordinates of r',x,y,z
         stop
      endif

      open(21, file='Exdx.mat', status='old', form='unformatted', access
     $     ='direct', recl=16)   
      open(22, file='Eydx.mat', status='old', form='unformatted', access
     $     ='direct', recl=16)
      open(23, file='Ezdx.mat', status='old', form='unformatted', access
     $     ='direct', recl=16)
      open(24, file='Exdy.mat', status='old', form='unformatted', access
     $     ='direct', recl=16)   
      open(25, file='Eydy.mat', status='old', form='unformatted', access
     $     ='direct', recl=16)
      open(26, file='Ezdy.mat', status='old', form='unformatted', access
     $     ='direct', recl=16)
      open(27, file='Exdz.mat', status='old', form='unformatted', access
     $     ='direct', recl=16)   
      open(28, file='Eydz.mat', status='old', form='unformatted', access
     $     ='direct', recl=16)
      open(29, file='Ezdz.mat', status='old', form='unformatted', access
     $     ='direct', recl=16)


      imin=IDNINT(dint((x-xmin)/dx))+1
      jmin=IDNINT(dint((y-ymin)/dy))+1
      kmin=IDNINT(dint((z-zmin)/dz))+1
c     write(*,*) 'ijk derive',imin,jmin,kmin
      xd=(x-xmin)/dx-dble(imin-1)
      yd=(y-ymin)/dy-dble(jmin-1)
      zd=(z-zmin)/dz-dble(kmin-1)
      

      if (test.eq.4) then

         if (xd.le.dx/1000.d0.and.yd.le.dy/1000.d0.and.zd.le.dz/1000.d0)
     $        then
c     on tombe sur une maille
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Ederint(1,1)
            read(22,rec=ii) Ederint(2,1)
            read(23,rec=ii) Ederint(3,1)
            read(24,rec=ii) Ederint(1,2)
            read(25,rec=ii) Ederint(2,2)
            read(26,rec=ii) Ederint(3,2)
            read(27,rec=ii) Ederint(1,3)
            read(28,rec=ii) Ederint(2,3)
            read(29,rec=ii) Ederint(3,3)
c     write(*,*) 'ii',ii,'Eder',Ederint
         elseif (xd.le.dx/1000.d0.and.yd.le.dy/1000.d0) then

            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(1,1,1)
            read(22,rec=ii) Eder(1,2,1)
            read(23,rec=ii) Eder(1,3,1)
            read(24,rec=ii) Eder(1,1,2)
            read(25,rec=ii) Eder(1,2,2)
            read(26,rec=ii) Eder(1,3,2)
            read(27,rec=ii) Eder(1,1,3)
            read(28,rec=ii) Eder(1,2,3)
            read(29,rec=ii) Eder(1,3,3)
            
            ii=imin+nx*(jmin-1)+nxy*kmin
            read(21,rec=ii) Eder(5,1,1)
            read(22,rec=ii) Eder(5,2,1)
            read(23,rec=ii) Eder(5,3,1)
            read(24,rec=ii) Eder(5,1,2)
            read(25,rec=ii) Eder(5,2,2)
            read(26,rec=ii) Eder(5,3,2)
            read(27,rec=ii) Eder(5,1,3)
            read(28,rec=ii) Eder(5,2,3)
            read(29,rec=ii) Eder(5,3,3)

            do kk=1,3
               do jj=1,3
                  Ederint(kk,jj)=(1.d0-zd)*Eder(1,kk,jj)+zd*Eder(5,kk
     $                 ,jj)
               enddo
            enddo
         elseif (xd.le.dx/1000.d0.and.zd.le.dz/1000.d0) then
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(1,1,1)
            read(22,rec=ii) Eder(1,2,1)
            read(23,rec=ii) Eder(1,3,1)
            read(24,rec=ii) Eder(1,1,2)
            read(25,rec=ii) Eder(1,2,2)
            read(26,rec=ii) Eder(1,3,2)
            read(27,rec=ii) Eder(1,1,3)
            read(28,rec=ii) Eder(1,2,3)
            read(29,rec=ii) Eder(1,3,3)

            ii=imin+nx*jmin+nxy*(kmin-1)
            read(21,rec=ii) Eder(3,1,1)
            read(22,rec=ii) Eder(3,2,1)
            read(23,rec=ii) Eder(3,3,1)
            read(24,rec=ii) Eder(3,1,2)
            read(25,rec=ii) Eder(3,2,2)
            read(26,rec=ii) Eder(3,3,2)
            read(27,rec=ii) Eder(3,1,3)
            read(28,rec=ii) Eder(3,2,3)
            read(29,rec=ii) Eder(3,3,3)

            do kk=1,3
               do jj=1,3
                  Ederint(kk,jj)=(1.d0-yd)*Eder(1,kk,jj)+yd*Eder(3,kk
     $                 ,jj)
               enddo
            enddo

         elseif (yd.le.dy/1000.d0.and.zd.le.dz/1000.d0) then
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(1,1,1)
            read(22,rec=ii) Eder(1,2,1)
            read(23,rec=ii) Eder(1,3,1)
            read(24,rec=ii) Eder(1,1,2)
            read(25,rec=ii) Eder(1,2,2)
            read(26,rec=ii) Eder(1,3,2)
            read(27,rec=ii) Eder(1,1,3)
            read(28,rec=ii) Eder(1,2,3)
            read(29,rec=ii) Eder(1,3,3)

            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(2,1,1)
            read(22,rec=ii) Eder(2,2,1)
            read(23,rec=ii) Eder(2,3,1)
            read(24,rec=ii) Eder(2,1,2)
            read(25,rec=ii) Eder(2,2,2)
            read(26,rec=ii) Eder(2,3,2)
            read(27,rec=ii) Eder(2,1,3)
            read(28,rec=ii) Eder(2,2,3)
            read(29,rec=ii) Eder(2,3,3)

            do kk=1,3
               do jj=1,3
                  Ederint(kk,jj)=(1.d0-xd)*Eder(1,kk,jj)+xd*Eder(2,kk
     $                 ,jj)
               enddo
            enddo

         elseif (xd.le.dx/1000.d0) then
c     1 coin
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(1,1,1)
            read(22,rec=ii) Eder(1,2,1)
            read(23,rec=ii) Eder(1,3,1)
            read(24,rec=ii) Eder(1,1,2)
            read(25,rec=ii) Eder(1,2,2)
            read(26,rec=ii) Eder(1,3,2)
            read(27,rec=ii) Eder(1,1,3)
            read(28,rec=ii) Eder(1,2,3)
            read(29,rec=ii) Eder(1,3,3)
c     3 coin
            ii=imin+nx*jmin+nxy*(kmin-1)
            read(21,rec=ii) Eder(3,1,1)
            read(22,rec=ii) Eder(3,2,1)
            read(23,rec=ii) Eder(3,3,1)
            read(24,rec=ii) Eder(3,1,2)
            read(25,rec=ii) Eder(3,2,2)
            read(26,rec=ii) Eder(3,3,2)
            read(27,rec=ii) Eder(3,1,3)
            read(28,rec=ii) Eder(3,2,3)
            read(29,rec=ii) Eder(3,3,3)
c     5 coin
            ii=imin+nx*(jmin-1)+nxy*kmin
            read(21,rec=ii) Eder(5,1,1)
            read(22,rec=ii) Eder(5,2,1)
            read(23,rec=ii) Eder(5,3,1)
            read(24,rec=ii) Eder(5,1,2)
            read(25,rec=ii) Eder(5,2,2)
            read(26,rec=ii) Eder(5,3,2)
            read(27,rec=ii) Eder(5,1,3)
            read(28,rec=ii) Eder(5,2,3)
            read(29,rec=ii) Eder(5,3,3)
c     7 coin
            ii=imin+nx*jmin+nxy*kmin
            read(21,rec=ii) Eder(7,1,1)
            read(22,rec=ii) Eder(7,2,1)
            read(23,rec=ii) Eder(7,3,1)
            read(24,rec=ii) Eder(7,1,2)
            read(25,rec=ii) Eder(7,2,2)
            read(26,rec=ii) Eder(7,3,2)
            read(27,rec=ii) Eder(7,1,3)
            read(28,rec=ii) Eder(7,2,3)
            read(29,rec=ii) Eder(7,3,3)

            do kk=1,3
               do jj=1,3
                  Ederint(kk,jj)=(1.d0-zd)*((1-yd)*Eder(1,kk,jj)+yd
     $                 *Eder(3,kk,jj))+zd*(Eder(5,kk,jj)*(1.d0-yd)
     $                 +Eder(7,kk,jj)*yd)
               enddo
            enddo

         elseif (yd.le.dy/1000.d0) then

            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(1,1,1)
            read(22,rec=ii) Eder(1,2,1)
            read(23,rec=ii) Eder(1,3,1)
            read(24,rec=ii) Eder(1,1,2)
            read(25,rec=ii) Eder(1,2,2)
            read(26,rec=ii) Eder(1,3,2)
            read(27,rec=ii) Eder(1,1,3)
            read(28,rec=ii) Eder(1,2,3)
            read(29,rec=ii) Eder(1,3,3)

            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(2,1,1)
            read(22,rec=ii) Eder(2,2,1)
            read(23,rec=ii) Eder(2,3,1)
            read(24,rec=ii) Eder(2,1,2)
            read(25,rec=ii) Eder(2,2,2)
            read(26,rec=ii) Eder(2,3,2)
            read(27,rec=ii) Eder(2,1,3)
            read(28,rec=ii) Eder(2,2,3)
            read(29,rec=ii) Eder(2,3,3)

            ii=imin+nx*(jmin-1)+nxy*kmin
            read(21,rec=ii) Eder(5,1,1)
            read(22,rec=ii) Eder(5,2,1)
            read(23,rec=ii) Eder(5,3,1)
            read(24,rec=ii) Eder(5,1,2)
            read(25,rec=ii) Eder(5,2,2)
            read(26,rec=ii) Eder(5,3,2)
            read(27,rec=ii) Eder(5,1,3)
            read(28,rec=ii) Eder(5,2,3)
            read(29,rec=ii) Eder(5,3,3)

            ii=imin+1+nx*(jmin-1)+nxy*kmin
            read(21,rec=ii) Eder(6,1,1)
            read(22,rec=ii) Eder(6,2,1)
            read(23,rec=ii) Eder(6,3,1)
            read(24,rec=ii) Eder(6,1,2)
            read(25,rec=ii) Eder(6,2,2)
            read(26,rec=ii) Eder(6,3,2)
            read(27,rec=ii) Eder(6,1,3)
            read(28,rec=ii) Eder(6,2,3)
            read(29,rec=ii) Eder(6,3,3)


            do kk=1,3
               do jj=1,3
                  Ederint(kk,jj)=(1.d0-zd)*((Eder(1,kk,jj)*(1.d0-xd)
     $                 +Eder(2,kk,jj)*xd))+zd*((Eder(5,kk,jj)*(1.d0-xd)
     $                 +Eder(6,kk,jj)*xd))
               enddo
            enddo
         elseif (zd.le.dz/1000.d0) then
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(1,1,1)
            read(22,rec=ii) Eder(1,2,1)
            read(23,rec=ii) Eder(1,3,1)
            read(24,rec=ii) Eder(1,1,2)
            read(25,rec=ii) Eder(1,2,2)
            read(26,rec=ii) Eder(1,3,2)
            read(27,rec=ii) Eder(1,1,3)
            read(28,rec=ii) Eder(1,2,3)
            read(29,rec=ii) Eder(1,3,3)

            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(2,1,1)
            read(22,rec=ii) Eder(2,2,1)
            read(23,rec=ii) Eder(2,3,1)
            read(24,rec=ii) Eder(2,1,2)
            read(25,rec=ii) Eder(2,2,2)
            read(26,rec=ii) Eder(2,3,2)
            read(27,rec=ii) Eder(2,1,3)
            read(28,rec=ii) Eder(2,2,3)
            read(29,rec=ii) Eder(2,3,3)

            ii=imin+nx*jmin+nxy*(kmin-1)
            read(21,rec=ii) Eder(3,1,1)
            read(22,rec=ii) Eder(3,2,1)
            read(23,rec=ii) Eder(3,3,1)
            read(24,rec=ii) Eder(3,1,2)
            read(25,rec=ii) Eder(3,2,2)
            read(26,rec=ii) Eder(3,3,2)
            read(27,rec=ii) Eder(3,1,3)
            read(28,rec=ii) Eder(3,2,3)
            read(29,rec=ii) Eder(3,3,3)

            ii=imin+1+nx*jmin+nxy*(kmin-1)
            read(21,rec=ii) Eder(4,1,1)
            read(22,rec=ii) Eder(4,2,1)
            read(23,rec=ii) Eder(4,3,1)
            read(24,rec=ii) Eder(4,1,2)
            read(25,rec=ii) Eder(4,2,2)
            read(26,rec=ii) Eder(4,3,2)
            read(27,rec=ii) Eder(4,1,3)
            read(28,rec=ii) Eder(4,2,3)
            read(29,rec=ii) Eder(4,3,3)

            do kk=1,3
               do jj=1,3
                  Ederint(kk,jj)=(1.d0-yd)*(Eder(1,kk,jj)*(1.d0-xd)
     $                 +Eder(2,kk,jj)*xd)+yd*(Eder(3,kk,jj)*(1.d0-xd) +
     $                 Eder(4,kk,jj)*xd)
               enddo
            enddo
         else

            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(1,1,1)
            read(22,rec=ii) Eder(1,2,1)
            read(23,rec=ii) Eder(1,3,1)
            read(24,rec=ii) Eder(1,1,2)
            read(25,rec=ii) Eder(1,2,2)
            read(26,rec=ii) Eder(1,3,2)
            read(27,rec=ii) Eder(1,1,3)
            read(28,rec=ii) Eder(1,2,3)
            read(29,rec=ii) Eder(1,3,3)

            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(2,1,1)
            read(22,rec=ii) Eder(2,2,1)
            read(23,rec=ii) Eder(2,3,1)
            read(24,rec=ii) Eder(2,1,2)
            read(25,rec=ii) Eder(2,2,2)
            read(26,rec=ii) Eder(2,3,2)
            read(27,rec=ii) Eder(2,1,3)
            read(28,rec=ii) Eder(2,2,3)
            read(29,rec=ii) Eder(2,3,3)

            ii=imin+nx*jmin+nxy*(kmin-1)
            read(21,rec=ii) Eder(3,1,1)
            read(22,rec=ii) Eder(3,2,1)
            read(23,rec=ii) Eder(3,3,1)
            read(24,rec=ii) Eder(3,1,2)
            read(25,rec=ii) Eder(3,2,2)
            read(26,rec=ii) Eder(3,3,2)
            read(27,rec=ii) Eder(3,1,3)
            read(28,rec=ii) Eder(3,2,3)
            read(29,rec=ii) Eder(3,3,3)

            ii=imin+1+nx*jmin+nxy*(kmin-1)
            read(21,rec=ii) Eder(4,1,1)
            read(22,rec=ii) Eder(4,2,1)
            read(23,rec=ii) Eder(4,3,1)
            read(24,rec=ii) Eder(4,1,2)
            read(25,rec=ii) Eder(4,2,2)
            read(26,rec=ii) Eder(4,3,2)
            read(27,rec=ii) Eder(4,1,3)
            read(28,rec=ii) Eder(4,2,3)
            read(29,rec=ii) Eder(4,3,3)

            ii=imin+nx*(jmin-1)+nxy*kmin
            read(21,rec=ii) Eder(5,1,1)
            read(22,rec=ii) Eder(5,2,1)
            read(23,rec=ii) Eder(5,3,1)
            read(24,rec=ii) Eder(5,1,2)
            read(25,rec=ii) Eder(5,2,2)
            read(26,rec=ii) Eder(5,3,2)
            read(27,rec=ii) Eder(5,1,3)
            read(28,rec=ii) Eder(5,2,3)
            read(29,rec=ii) Eder(5,3,3)

            ii=imin+1+nx*(jmin-1)+nxy*kmin
            read(21,rec=ii) Eder(6,1,1)
            read(22,rec=ii) Eder(6,2,1)
            read(23,rec=ii) Eder(6,3,1)
            read(24,rec=ii) Eder(6,1,2)
            read(25,rec=ii) Eder(6,2,2)
            read(26,rec=ii) Eder(6,3,2)
            read(27,rec=ii) Eder(6,1,3)
            read(28,rec=ii) Eder(6,2,3)
            read(29,rec=ii) Eder(6,3,3)

            ii=imin+nx*jmin+nxy*kmin
            read(21,rec=ii) Eder(7,1,1)
            read(22,rec=ii) Eder(7,2,1)
            read(23,rec=ii) Eder(7,3,1)
            read(24,rec=ii) Eder(7,1,2)
            read(25,rec=ii) Eder(7,2,2)
            read(26,rec=ii) Eder(7,3,2)
            read(27,rec=ii) Eder(7,1,3)
            read(28,rec=ii) Eder(7,2,3)
            read(29,rec=ii) Eder(7,3,3)

            ii=imin+1+nx*jmin+nxy*kmin
            read(21,rec=ii) Eder(8,1,1)
            read(22,rec=ii) Eder(8,2,1)
            read(23,rec=ii) Eder(8,3,1)
            read(24,rec=ii) Eder(8,1,2)
            read(25,rec=ii) Eder(8,2,2)
            read(26,rec=ii) Eder(8,3,2)
            read(27,rec=ii) Eder(8,1,3)
            read(28,rec=ii) Eder(8,2,3)
            read(29,rec=ii) Eder(8,3,3)

c     interpolation trilineaire      
            do kk=1,3
               do jj=1,3
                  Ederint(kk,jj)=(1.d0-zd)*((1-yd)*(Eder(1,kk,jj)*(1.d0
     $                 -xd)+Eder(2,kk,jj)*xd)+yd*(Eder(3,kk,jj)*(1.d0
     $                 -xd)+Eder(4,kk,jj)*xd))+zd*((Eder(5,kk,jj)*(1.d0
     $                 -xd)+Eder(6,kk,jj)*xd)*(1.d0-yd)+(Eder(7,kk,jj)
     $                 *(1.d0-xd)+Eder(8,kk,jj)*xd)*yd)
               enddo
            enddo
         endif

      elseif (test.eq.1) then
         if (xd.le.dx/1000.d0.and.yd.le.dy/1000.d0.and.zd.le.dz/1000.d0)
     $        then
c     on tombe sur une maille
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Ederint(1,1)
            read(22,rec=ii) Ederint(2,1)
            read(23,rec=ii) Ederint(3,1)

         elseif (xd.le.dx/1000.d0.and.yd.le.dy/1000.d0) then

            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(1,1,1)
            read(22,rec=ii) Eder(1,2,1)
            read(23,rec=ii) Eder(1,3,1)
            ii=imin+nx*(jmin-1)+nxy*kmin
            read(21,rec=ii) Eder(5,1,1)
            read(22,rec=ii) Eder(5,2,1)
            read(23,rec=ii) Eder(5,3,1)
            do kk=1,3
               Ederint(kk,1)=(1.d0-zd)*Eder(1,kk,1)+zd*Eder(5,kk,1)
            enddo

         elseif (xd.le.dx/1000.d0.and.zd.le.dz/1000.d0) then
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(1,1,1)
            read(22,rec=ii) Eder(1,2,1)
            read(23,rec=ii) Eder(1,3,1)
            ii=imin+nx*jmin+nxy*(kmin-1)
            read(21,rec=ii) Eder(3,1,1)
            read(22,rec=ii) Eder(3,2,1)
            read(23,rec=ii) Eder(3,3,1)
        
            do kk=1,3             
               Ederint(kk,1)=(1.d0-yd)*Eder(1,kk,1)+yd*Eder(3,kk,1)
            enddo

         elseif (yd.le.dy/1000.d0.and.zd.le.dz/1000.d0) then
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(1,1,1)
            read(22,rec=ii) Eder(1,2,1)
            read(23,rec=ii) Eder(1,3,1)
            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(2,1,1)
            read(22,rec=ii) Eder(2,2,1)
            read(23,rec=ii) Eder(2,3,1)

            do kk=1,3
               Ederint(kk,1)=(1.d0-xd)*Eder(1,kk,1)+xd*Eder(2,kk,1)
            enddo

         elseif (xd.le.dx/1000.d0) then
c     1 coin
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(1,1,1)
            read(22,rec=ii) Eder(1,2,1)
            read(23,rec=ii) Eder(1,3,1)
c     3 coin
            ii=imin+nx*jmin+nxy*(kmin-1)
            read(21,rec=ii) Eder(3,1,1)
            read(22,rec=ii) Eder(3,2,1)
            read(23,rec=ii) Eder(3,3,1)
c     5 coin
            ii=imin+nx*(jmin-1)+nxy*kmin
            read(21,rec=ii) Eder(5,1,1)
            read(22,rec=ii) Eder(5,2,1)
            read(23,rec=ii) Eder(5,3,1)
c     7 coin
            ii=imin+nx*jmin+nxy*kmin
            read(21,rec=ii) Eder(7,1,1)
            read(22,rec=ii) Eder(7,2,1)
            read(23,rec=ii) Eder(7,3,1)

            do kk=1,3
               Ederint(kk,1)=(1.d0-zd)*((1-yd)*Eder(1,kk,1)+yd*Eder(3
     $              ,kk,1))+zd*(Eder(5,kk,1)*(1.d0-yd)+Eder(7,kk,1)*yd)
            enddo

         elseif (yd.le.dy/1000.d0) then

            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(1,1,1)
            read(22,rec=ii) Eder(1,2,1)
            read(23,rec=ii) Eder(1,3,1)
            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(2,1,1)
            read(22,rec=ii) Eder(2,2,1)
            read(23,rec=ii) Eder(2,3,1)          
            ii=imin+nx*(jmin-1)+nxy*kmin
            read(21,rec=ii) Eder(5,1,1)
            read(22,rec=ii) Eder(5,2,1)
            read(23,rec=ii) Eder(5,3,1)
            ii=imin+1+nx*(jmin-1)+nxy*kmin
            read(21,rec=ii) Eder(6,1,1)
            read(22,rec=ii) Eder(6,2,1)
            read(23,rec=ii) Eder(6,3,1)

            do kk=1,3
               Ederint(kk,1)=(1.d0-zd)*((Eder(1,kk,1)*(1.d0-xd)+Eder(2
     $              ,kk,1)*xd))+zd*((Eder(5,kk,1)*(1.d0-xd)+Eder(6,kk
     $              ,1)*xd))
            enddo
         elseif (zd.le.dz/1000.d0) then
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(1,1,1)
            read(22,rec=ii) Eder(1,2,1)
            read(23,rec=ii) Eder(1,3,1)
            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(2,1,1)
            read(22,rec=ii) Eder(2,2,1)
            read(23,rec=ii) Eder(2,3,1)
            ii=imin+nx*jmin+nxy*(kmin-1)
            read(21,rec=ii) Eder(3,1,1)
            read(22,rec=ii) Eder(3,2,1)
            read(23,rec=ii) Eder(3,3,1)
            ii=imin+1+nx*jmin+nxy*(kmin-1)
            read(21,rec=ii) Eder(4,1,1)
            read(22,rec=ii) Eder(4,2,1)
            read(23,rec=ii) Eder(4,3,1)
            do kk=1,3
               Ederint(kk,1)=(1.d0-yd)*(Eder(1,kk,1)*(1.d0-xd)+Eder(2
     $              ,kk,1)*xd)+yd*(Eder(3,kk,1)*(1.d0-xd)+Eder(4,kk,1)
     $              *xd)
            enddo
         else

            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(1,1,1)
            read(22,rec=ii) Eder(1,2,1)
            read(23,rec=ii) Eder(1,3,1)
            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(21,rec=ii) Eder(2,1,1)
            read(22,rec=ii) Eder(2,2,1)
            read(23,rec=ii) Eder(2,3,1)
            ii=imin+nx*jmin+nxy*(kmin-1)
            read(21,rec=ii) Eder(3,1,1)
            read(22,rec=ii) Eder(3,2,1)
            read(23,rec=ii) Eder(3,3,1)
            ii=imin+1+nx*jmin+nxy*(kmin-1)
            read(21,rec=ii) Eder(4,1,1)
            read(22,rec=ii) Eder(4,2,1)
            read(23,rec=ii) Eder(4,3,1)
            ii=imin+nx*(jmin-1)+nxy*kmin
            read(21,rec=ii) Eder(5,1,1)
            read(22,rec=ii) Eder(5,2,1)
            read(23,rec=ii) Eder(5,3,1)
            ii=imin+1+nx*(jmin-1)+nxy*kmin
            read(21,rec=ii) Eder(6,1,1)
            read(22,rec=ii) Eder(6,2,1)
            read(23,rec=ii) Eder(6,3,1)
            ii=imin+nx*jmin+nxy*kmin
            read(21,rec=ii) Eder(7,1,1)
            read(22,rec=ii) Eder(7,2,1)
            read(23,rec=ii) Eder(7,3,1)
            ii=imin+1+nx*jmin+nxy*kmin
            read(21,rec=ii) Eder(8,1,1)
            read(22,rec=ii) Eder(8,2,1)
            read(23,rec=ii) Eder(8,3,1)

c     interpolation trilineaire      
            do kk=1,3            
               Ederint(kk,1)=(1.d0-zd)*((1-yd)*(Eder(1,kk,1)*(1.d0-xd)
     $              +Eder(2,kk,1)*xd)+yd*(Eder(3,kk,1)*(1.d0-xd)+Eder(4
     $              ,kk,1)*xd))+zd*((Eder(5,kk,1)*(1.d0-xd)+Eder(6,kk
     $              ,1)*xd)*(1.d0-yd)+(Eder(7,kk,1)*(1.d0-xd)+Eder(8,kk
     $              ,1)*xd)*yd)
            enddo
         endif

      elseif (test.eq.2) then
         if (xd.le.dx/1000.d0.and.yd.le.dy/1000.d0.and.zd.le.dz/1000.d0)
     $        then
c     on tombe sur une maille
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)           
            read(24,rec=ii) Ederint(1,2)
            read(25,rec=ii) Ederint(2,2)
            read(26,rec=ii) Ederint(3,2)           
         elseif (xd.le.dx/1000.d0.and.yd.le.dy/1000.d0) then

            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(24,rec=ii) Eder(1,1,2)
            read(25,rec=ii) Eder(1,2,2)
            read(26,rec=ii) Eder(1,3,2)
            ii=imin+nx*(jmin-1)+nxy*kmin
            read(24,rec=ii) Eder(5,1,2)
            read(25,rec=ii) Eder(5,2,2)
            read(26,rec=ii) Eder(5,3,2)

            do kk=1,3
               Ederint(kk,2)=(1.d0-zd)*Eder(1,kk,2)+zd*Eder(5,kk,2)
            enddo
         elseif (xd.le.dx/1000.d0.and.zd.le.dz/1000.d0) then
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(24,rec=ii) Eder(1,1,2)
            read(25,rec=ii) Eder(1,2,2)
            read(26,rec=ii) Eder(1,3,2)
            ii=imin+nx*jmin+nxy*(kmin-1)
            read(24,rec=ii) Eder(3,1,2)
            read(25,rec=ii) Eder(3,2,2)
            read(26,rec=ii) Eder(3,3,2)

            do kk=1,3
               Ederint(kk,2)=(1.d0-yd)*Eder(1,kk,2)+yd*Eder(3,kk,2)
            enddo

         elseif (yd.le.dy/1000.d0.and.zd.le.dz/1000.d0) then
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(24,rec=ii) Eder(1,1,2)
            read(25,rec=ii) Eder(1,2,2)
            read(26,rec=ii) Eder(1,3,2)
            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(24,rec=ii) Eder(2,1,2)
            read(25,rec=ii) Eder(2,2,2)
            read(26,rec=ii) Eder(2,3,2)

            do kk=1,3
               Ederint(kk,2)=(1.d0-xd)*Eder(1,kk,2)+xd*Eder(2,kk,2)
            enddo
         elseif (xd.le.dx/1000.d0) then
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(24,rec=ii) Eder(1,1,2)
            read(25,rec=ii) Eder(1,2,2)
            read(26,rec=ii) Eder(1,3,2)
            ii=imin+nx*jmin+nxy*(kmin-1)
            read(24,rec=ii) Eder(3,1,2)
            read(25,rec=ii) Eder(3,2,2)
            read(26,rec=ii) Eder(3,3,2)
            ii=imin+nx*(jmin-1)+nxy*kmin
            read(24,rec=ii) Eder(5,1,2)
            read(25,rec=ii) Eder(5,2,2)
            read(26,rec=ii) Eder(5,3,2)
            ii=imin+nx*jmin+nxy*kmin
            read(24,rec=ii) Eder(7,1,2)
            read(25,rec=ii) Eder(7,2,2)
            read(26,rec=ii) Eder(7,3,2)

            do kk=1,3
               Ederint(kk,2)=(1.d0-zd)*((1-yd)*Eder(1,kk,2)+yd*Eder(3
     $              ,kk,2))+zd*(Eder(5,kk,2)*(1.d0-yd)+Eder(7,kk,2)*yd)
            enddo

         elseif (yd.le.dy/1000.d0) then

            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(24,rec=ii) Eder(1,1,2)
            read(25,rec=ii) Eder(1,2,2)
            read(26,rec=ii) Eder(1,3,2)
            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(24,rec=ii) Eder(2,1,2)
            read(25,rec=ii) Eder(2,2,2)
            read(26,rec=ii) Eder(2,3,2)
            ii=imin+nx*(jmin-1)+nxy*kmin
            read(24,rec=ii) Eder(5,1,2)
            read(25,rec=ii) Eder(5,2,2)
            read(26,rec=ii) Eder(5,3,2)
            ii=imin+1+nx*(jmin-1)+nxy*kmin
            read(24,rec=ii) Eder(6,1,2)
            read(25,rec=ii) Eder(6,2,2)
            read(26,rec=ii) Eder(6,3,2)

            do kk=1,3
               Ederint(kk,2)=(1.d0-zd)*((Eder(1,kk,2)*(1.d0-xd) +Eder(2
     $              ,kk,2)*xd))+zd*((Eder(5,kk,2)*(1.d0-xd)+Eder(6,kk
     $              ,2)*xd))
            enddo
         elseif (zd.le.dz/1000.d0) then
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(24,rec=ii) Eder(1,1,2)
            read(25,rec=ii) Eder(1,2,2)
            read(26,rec=ii) Eder(1,3,2)
            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(24,rec=ii) Eder(2,1,2)
            read(25,rec=ii) Eder(2,2,2)
            read(26,rec=ii) Eder(2,3,2)
            ii=imin+nx*jmin+nxy*(kmin-1)
            read(24,rec=ii) Eder(3,1,2)
            read(25,rec=ii) Eder(3,2,2)
            read(26,rec=ii) Eder(3,3,2)
            ii=imin+1+nx*jmin+nxy*(kmin-1)
            read(24,rec=ii) Eder(4,1,2)
            read(25,rec=ii) Eder(4,2,2)
            read(26,rec=ii) Eder(4,3,2)

            do kk=1,3
               Ederint(kk,2)=(1.d0-yd)*(Eder(1,kk,2)*(1.d0-xd)+Eder(2
     $              ,kk,2)*xd)+yd*(Eder(3,kk,2)*(1.d0-xd)+Eder(4,kk,2)
     $              *xd)
            enddo
         else

            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(24,rec=ii) Eder(1,1,2)
            read(25,rec=ii) Eder(1,2,2)
            read(26,rec=ii) Eder(1,3,2)
            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(24,rec=ii) Eder(2,1,2)
            read(25,rec=ii) Eder(2,2,2)
            read(26,rec=ii) Eder(2,3,2)
            ii=imin+nx*jmin+nxy*(kmin-1)
            read(24,rec=ii) Eder(3,1,2)
            read(25,rec=ii) Eder(3,2,2)
            read(26,rec=ii) Eder(3,3,2)
            ii=imin+1+nx*jmin+nxy*(kmin-1)
            read(24,rec=ii) Eder(4,1,2)
            read(25,rec=ii) Eder(4,2,2)
            read(26,rec=ii) Eder(4,3,2)
            ii=imin+nx*(jmin-1)+nxy*kmin
            read(24,rec=ii) Eder(5,1,2)
            read(25,rec=ii) Eder(5,2,2)
            read(26,rec=ii) Eder(5,3,2)
            ii=imin+1+nx*(jmin-1)+nxy*kmin
            read(24,rec=ii) Eder(6,1,2)
            read(25,rec=ii) Eder(6,2,2)
            read(26,rec=ii) Eder(6,3,2)
            ii=imin+nx*jmin+nxy*kmin
            read(24,rec=ii) Eder(7,1,2)
            read(25,rec=ii) Eder(7,2,2)
            read(26,rec=ii) Eder(7,3,2)
            ii=imin+1+nx*jmin+nxy*kmin
            read(24,rec=ii) Eder(8,1,2)
            read(25,rec=ii) Eder(8,2,2)
            read(26,rec=ii) Eder(8,3,2)

c     interpolation trilineaire      
            do kk=1,3
               Ederint(kk,2)=(1.d0-zd)*((1-yd)*(Eder(1,kk,2)*(1.d0-xd)
     $              +Eder(2,kk,2)*xd)+yd*(Eder(3,kk,2)*(1.d0-xd)+Eder(4
     $              ,kk,2)*xd))+zd*((Eder(5,kk,2)*(1.d0-xd)+Eder(6,kk
     $              ,2)*xd)*(1.d0-yd)+(Eder(7,kk,2)*(1.d0-xd)+Eder(8,kk
     $              ,2)*xd)*yd)
            enddo
         endif

      elseif (test.eq.3) then
         if (xd.le.dx/1000.d0.and.yd.le.dy/1000.d0.and.zd.le.dz/1000.d0)
     $        then
c     on tombe sur une maille
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(27,rec=ii) Ederint(1,3)
            read(28,rec=ii) Ederint(2,3)
            read(29,rec=ii) Ederint(3,3)

         elseif (xd.le.dx/1000.d0.and.yd.le.dy/1000.d0) then

            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(27,rec=ii) Eder(1,1,3)
            read(28,rec=ii) Eder(1,2,3)
            read(29,rec=ii) Eder(1,3,3)
            ii=imin+nx*(jmin-1)+nxy*kmin
            read(27,rec=ii) Eder(5,1,3)
            read(28,rec=ii) Eder(5,2,3)
            read(29,rec=ii) Eder(5,3,3)

            do kk=1,3
               Ederint(kk,3)=(1.d0-zd)*Eder(1,kk,3)+zd*Eder(5,kk,3)
            enddo

         elseif (xd.le.dx/1000.d0.and.zd.le.dz/1000.d0) then
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(27,rec=ii) Eder(1,1,3)
            read(28,rec=ii) Eder(1,2,3)
            read(29,rec=ii) Eder(1,3,3)
            ii=imin+nx*jmin+nxy*(kmin-1)
            read(27,rec=ii) Eder(3,1,3)
            read(28,rec=ii) Eder(3,2,3)
            read(29,rec=ii) Eder(3,3,3)

            do kk=1,3
               Ederint(kk,3)=(1.d0-yd)*Eder(1,kk,3)+yd*Eder(3,kk,3)
            enddo

         elseif (yd.le.dy/1000.d0.and.zd.le.dz/1000.d0) then
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(27,rec=ii) Eder(1,1,3)
            read(28,rec=ii) Eder(1,2,3)
            read(29,rec=ii) Eder(1,3,3)
            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(27,rec=ii) Eder(2,1,3)
            read(28,rec=ii) Eder(2,2,3)
            read(29,rec=ii) Eder(2,3,3)

            do kk=1,3
               Ederint(kk,3)=(1.d0-xd)*Eder(1,kk,3)+xd*Eder(2,kk,3)
            enddo

         elseif (xd.le.dx/1000.d0) then

            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(27,rec=ii) Eder(1,1,3)
            read(28,rec=ii) Eder(1,2,3)
            read(29,rec=ii) Eder(1,3,3)
            ii=imin+nx*jmin+nxy*(kmin-1)
            read(27,rec=ii) Eder(3,1,3)
            read(28,rec=ii) Eder(3,2,3)
            read(29,rec=ii) Eder(3,3,3)
            ii=imin+nx*(jmin-1)+nxy*kmin
            read(27,rec=ii) Eder(5,1,3)
            read(28,rec=ii) Eder(5,2,3)
            read(29,rec=ii) Eder(5,3,3)
            ii=imin+nx*jmin+nxy*kmin
            read(27,rec=ii) Eder(7,1,3)
            read(28,rec=ii) Eder(7,2,3)
            read(29,rec=ii) Eder(7,3,3)

            do kk=1,3
               Ederint(kk,3)=(1.d0-zd)*((1-yd)*Eder(1,kk,3)+yd*Eder(3
     $              ,kk,3))+zd*(Eder(5,kk,3)*(1.d0-yd)+Eder(7,kk,3)*yd)
            enddo

         elseif (yd.le.dy/1000.d0) then

            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(27,rec=ii) Eder(1,1,3)
            read(28,rec=ii) Eder(1,2,3)
            read(29,rec=ii) Eder(1,3,3)
            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(27,rec=ii) Eder(2,1,3)
            read(28,rec=ii) Eder(2,2,3)
            read(29,rec=ii) Eder(2,3,3)
            ii=imin+nx*(jmin-1)+nxy*kmin
            read(27,rec=ii) Eder(5,1,3)
            read(28,rec=ii) Eder(5,2,3)
            read(29,rec=ii) Eder(5,3,3)
            ii=imin+1+nx*(jmin-1)+nxy*kmin
            read(27,rec=ii) Eder(6,1,3)
            read(28,rec=ii) Eder(6,2,3)
            read(29,rec=ii) Eder(6,3,3)

            do kk=1,3
               Ederint(kk,3)=(1.d0-zd)*((Eder(1,kk,3)*(1.d0-xd)+Eder(2
     $              ,kk,3)*xd))+zd*((Eder(5,kk,3)*(1.d0-xd)+Eder(6,kk
     $              ,3)*xd))
            enddo

         elseif (zd.le.dz/1000.d0) then
            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(27,rec=ii) Eder(1,1,3)
            read(28,rec=ii) Eder(1,2,3)
            read(29,rec=ii) Eder(1,3,3)
            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(27,rec=ii) Eder(2,1,3)
            read(28,rec=ii) Eder(2,2,3)
            read(29,rec=ii) Eder(2,3,3)
            ii=imin+nx*jmin+nxy*(kmin-1)
            read(27,rec=ii) Eder(3,1,3)
            read(28,rec=ii) Eder(3,2,3)
            read(29,rec=ii) Eder(3,3,3)
            ii=imin+1+nx*jmin+nxy*(kmin-1)
            read(27,rec=ii) Eder(4,1,3)
            read(28,rec=ii) Eder(4,2,3)
            read(29,rec=ii) Eder(4,3,3)

            do kk=1,3
               Ederint(kk,3)=(1.d0-yd)*(Eder(1,kk,3)*(1.d0-xd)+Eder(2
     $              ,kk,3)*xd)+yd*(Eder(3,kk,3)*(1.d0-xd)+Eder(4,kk,3)
     $              *xd)
            enddo
         else

            ii=imin+nx*(jmin-1)+nxy*(kmin-1)
            read(27,rec=ii) Eder(1,1,3)
            read(28,rec=ii) Eder(1,2,3)
            read(29,rec=ii) Eder(1,3,3)
            ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
            read(27,rec=ii) Eder(2,1,3)
            read(28,rec=ii) Eder(2,2,3)
            read(29,rec=ii) Eder(2,3,3)
            ii=imin+nx*jmin+nxy*(kmin-1)
            read(27,rec=ii) Eder(3,1,3)
            read(28,rec=ii) Eder(3,2,3)
            read(29,rec=ii) Eder(3,3,3)
            ii=imin+1+nx*jmin+nxy*(kmin-1)
            read(27,rec=ii) Eder(4,1,3)
            read(28,rec=ii) Eder(4,2,3)
            read(29,rec=ii) Eder(4,3,3)
            ii=imin+nx*(jmin-1)+nxy*kmin
            read(27,rec=ii) Eder(5,1,3)
            read(28,rec=ii) Eder(5,2,3)
            read(29,rec=ii) Eder(5,3,3)
            ii=imin+1+nx*(jmin-1)+nxy*kmin
            read(27,rec=ii) Eder(6,1,3)
            read(28,rec=ii) Eder(6,2,3)
            read(29,rec=ii) Eder(6,3,3)
            ii=imin+nx*jmin+nxy*kmin
            read(27,rec=ii) Eder(7,1,3)
            read(28,rec=ii) Eder(7,2,3)
            read(29,rec=ii) Eder(7,3,3)
            ii=imin+1+nx*jmin+nxy*kmin
            read(27,rec=ii) Eder(8,1,3)
            read(28,rec=ii) Eder(8,2,3)
            read(29,rec=ii) Eder(8,3,3)

c     interpolation trilineaire      
            do kk=1,3
               Ederint(kk,3)=(1.d0-zd)*((1-yd)*(Eder(1,kk,3)*(1.d0 -xd)
     $              +Eder(2,kk,3)*xd)+yd*(Eder(3,kk,3)*(1.d0-xd)+Eder(4
     $              ,kk,3)*xd))+zd*((Eder(5,kk,3)*(1.d0-xd)+Eder(6,kk
     $              ,3)*xd)*(1.d0-yd)+(Eder(7,kk,3)*(1.d0-xd)+Eder(8,kk
     $              ,3)*xd)*yd)
            enddo
         endif

      endif

      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)
      close(29)

      end
c     **********************************************************
c     **********************************************************
c     **********************************************************
      subroutine incidentarbitrarypos(xs,ys,zs,aretecube,Exx,Eyy,Ezz
     $     ,nstop,namefile,infostr)
      implicit none
      integer ii,nx,ny,nz,nxy,imin,jmin ,kmin,nstop ,ierror
      double precision dx,dy,dz,xmin,ymin,zmin,xmax,ymax,zmax,exr,exi
     $     ,eyr,eyi,ezr,ezi,xd,yd,zd,aretecube,dmax

      double complex Ex(8),Ey(8),Ez(8)

      double precision xs,ys,zs
      double complex Exx,Eyy,Ezz,uncomp,icomp
      character(64) infostr
      character(64) namefile

      uncomp=(1.d0,0.d0)
      icomp=(0.d0,1.d0)

c     read the input file
      open(15,file=namefile,status='old',iostat=ierror) 
      if (ierror.ne.0) then
         infostr='arbitrary incident: name file does not exist'
         nstop=1
         return
      endif

      read(15,*) nx,ny,nz           
      read(15,*) dx,dy,dz
      read(15,*) xmin,ymin,zmin
      close(15)
      nxy=nx*ny

      dx=dx*1.d-9
      dy=dy*1.d-9
      dz=dz*1.d-9
      dmax=max(dx,dy,dz)
      if (aretecube.le.dmax) then
         infostr='object mesh size < mesh size incident field'
         nstop=1
         return
      endif


      xmin=xmin*1.d-9
      ymin=ymin*1.d-9
      zmin=zmin*1.d-9

      xmax=xmin+dble(nx-1)*dx
      ymax=ymin+dble(ny-1)*dy
      zmax=zmin+dble(nz-1)*dz

c     verifier que l'objet est dans la boite
c      write(*,*) 'rou',xs,ys,zs
      if (xs.lt.xmin.or.xs.gt.xmax .or. ys.lt.ymin.or.ys.gt.ymax .or.
     $     zs.lt.zmin.or.zs.gt.zmax) then
         nstop=1
         infostr='near field not inside the incident box'
         write(*,*) 'objet not inside the incident box'
         write(*,*) xmin,xmax,ymin ,ymax,zmin,zmax
         write(*,*) xs,ys,zs
         return
      endif
c     ecriture des fichiers sous cette structure
c     do i=1,nz
c     do j=1,ny
c     do k=1,nx
c     x=x0+dble(k-1)*dx
c     y=y0+dble(j-1)*dy
c     z=z0+dble(i-1)*dz

      open(11, file='Exr.mat', status='old', form='formatted', access
     $     ='direct', recl=22,iostat=ierror) 
      if (ierror.ne.0) then
         infostr='arbitrary incident: Exr.mat does not exist'
         nstop=1
         return
      endif
      open(12, file='Exi.mat', status='old', form='formatted', access
     $     ='direct', recl=22,iostat=ierror)
      if (ierror.ne.0) then
         infostr='arbitrary incident: Exi.mat does not exist'
         nstop=1
         return
      endif
      open(13, file='Eyr.mat', status='old', form='formatted', access
     $     ='direct', recl=22,iostat=ierror)
      if (ierror.ne.0) then
         infostr='arbitrary incident: Eyr.mat does not exist'
         nstop=1
         return
      endif
      open(14, file='Eyi.mat', status='old', form='formatted', access
     $     ='direct', recl=22,iostat=ierror)
      if (ierror.ne.0) then
         infostr='arbitrary incident: Eyi.mat does not exist'
         nstop=1
         return
      endif
      open(15, file='Ezr.mat', status='old', form='formatted', access
     $     ='direct', recl=22,iostat=ierror)
      if (ierror.ne.0) then
         infostr='arbitrary incident: Ezr.mat does not exist'
         nstop=1
         return
      endif
      open(16, file='Ezi.mat', status='old', form='formatted', access
     $     ='direct', recl=22,iostat=ierror)
      if (ierror.ne.0) then
         infostr='arbitrary incident: Ezi.mat does not exist'
         nstop=1
         return
      endif

      imin=IDNINT(dint((xs-xmin)/dx))+1
      jmin=IDNINT(dint((ys-ymin)/dy))+1
      kmin=IDNINT(dint((zs-zmin)/dz))+1
      
      xd=(xs-xmin)/dx-dble(imin-1)
      yd=(ys-ymin)/dy-dble(jmin-1)
      zd=(zs-zmin)/dz-dble(kmin-1)
      
      if (xd.le.dx/1000.d0.and.yd.le.dy/1000.d0.and.zd.le.dz/1000.d0)
     $     then
c     on tombe sur une maille
         ii=imin+nx*(jmin-1)+nxy*(kmin-1)
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(1)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(1)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(1)=Ezr*uncomp+icomp*Ezi            
         Exx=Ex(1)
         Eyy=Ey(1)
         Ezz=Ez(1)

      elseif (xd.le.dx/1000.d0.and.yd.le.dy/1000.d0) then

         ii=imin+nx*(jmin-1)+nxy*(kmin-1)
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(1)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(1)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(1)=Ezr*uncomp+icomp*Ezi   

         ii=imin+nx*(jmin-1)+nxy*kmin
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(5)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(5)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(5)=Ezr*uncomp+icomp*Ezi

         Exx=(1.d0-zd)*Ex(1)+zd*Ex(5)
         Eyy=(1.d0-zd)*Ey(1)+zd*Ey(5)
         Ezz=(1.d0-zd)*Ez(1)+zd*Ez(5)
         

      elseif (xd.le.dx/1000.d0.and.zd.le.dz/1000.d0) then
         ii=imin+nx*(jmin-1)+nxy*(kmin-1)
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(1)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(1)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(1)=Ezr*uncomp+icomp*Ezi   

         ii=imin+nx*jmin+nxy*(kmin-1)
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(3)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(3)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(3)=Ezr*uncomp+icomp*Ezi

         Exx=(1.d0-yd)*Ex(1)+yd*Ex(3)
         Eyy=(1.d0-yd)*Ey(1)+yd*Ey(3)
         Ezz=(1.d0-yd)*Ez(1)+yd*Ez(3)

         

      elseif (yd.le.dy/1000.d0.and.zd.le.dz/1000.d0) then
         ii=imin+nx*(jmin-1)+nxy*(kmin-1)
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(1)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(1)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(1)=Ezr*uncomp+icomp*Ezi   

         ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(2)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(2)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(2)=Ezr*uncomp+icomp*Ezi

         Exx=Ex(1)*(1.d0-xd)+Ex(2)*xd
         Eyy=Ey(1)*(1.d0-xd)+Ey(2)*xd
         Ezz=Ez(1)*(1.d0-xd)+Ez(2)*xd

      elseif (xd.le.dx/1000.d0) then


c     1 coin
         ii=imin+nx*(jmin-1)+nxy*(kmin-1)
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(1)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(1)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(1)=Ezr*uncomp+icomp*Ezi

c     3 coin
         ii=imin+nx*jmin+nxy*(kmin-1)
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(3)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(3)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(3)=Ezr*uncomp+icomp*Ezi

c     5 coin
         ii=imin+nx*(jmin-1)+nxy*kmin
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(5)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(5)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(5)=Ezr*uncomp+icomp*Ezi


c     7 coin
         ii=imin+nx*jmin+nxy*kmin
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(7)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(7)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(7)=Ezr*uncomp+icomp*Ezi


         Exx=(1.d0-zd)*((1-yd)*Ex(1)+yd*Ex(3))+zd*(Ex(5) *(1.d0-yd)
     $        +Ex(7)*yd)
         Eyy=(1.d0-zd)*((1-yd)*Ey(1)+yd*Ey(3))+zd*(Ey(5) *(1.d0-yd)
     $        +Ey(7)*yd)
         Ezz=(1.d0-zd)*((1-yd)*Ez(1)+yd*Ez(3))+zd*(Ez(5) *(1.d0-yd)
     $        +Ez(7)*yd)

         
      elseif (yd.le.dy/1000.d0) then

         ii=imin+nx*(jmin-1)+nxy*(kmin-1)
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(1)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(1)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(1)=Ezr*uncomp+icomp*Ezi

         ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(2)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(2)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(2)=Ezr*uncomp+icomp*Ezi

         ii=imin+nx*(jmin-1)+nxy*kmin
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(5)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(5)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(5)=Ezr*uncomp+icomp*Ezi

         ii=imin+1+nx*(jmin-1)+nxy*kmin
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(6)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(6)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(6)=Ezr*uncomp+icomp*Ezi

         Exx=(1.d0-zd)*((Ex(1)*(1.d0-xd)+Ex(2)*xd)) +zd *((Ex(5) *(1.d0
     $        -xd)+Ex(6)*xd))
         Eyy=(1.d0-zd)*((Ey(1)*(1.d0-xd)+Ey(2)*xd)) +zd *((Ey(5) *(1.d0
     $        -xd)+Ey(6)*xd))
         Ezz=(1.d0-zd)*((Ez(1)*(1.d0-xd)+Ez(2)*xd)) +zd *((Ez(5) *(1.d0
     $        -xd)+Ez(6)*xd))


      elseif (zd.le.dz/1000.d0) then
c     1 coin
         ii=imin+nx*(jmin-1)+nxy*(kmin-1)
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(1)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(1)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(1)=Ezr*uncomp+icomp*Ezi

c     2 coin
         ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(2)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(2)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(2)=Ezr*uncomp+icomp*Ezi

c     3 coin
         ii=imin+nx*jmin+nxy*(kmin-1)
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(3)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(3)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(3)=Ezr*uncomp+icomp*Ezi

c     4 coin
         ii=imin+1+nx*jmin+nxy*(kmin-1)
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(4)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(4)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(4)=Ezr*uncomp+icomp*Ezi

         Exx=(1.d0-yd)*(Ex(1)*(1.d0-xd)+Ex(2)*xd)+yd*(Ex(3) *(1.d0-xd) +
     $        Ex(4)*xd)
         Eyy=(1.d0-yd)*(Ey(1)*(1.d0-xd)+Ey(2)*xd)+yd*(Ey(3) *(1.d0-xd)
     $        +Ey(4)*xd)
         Ezz=(1.d0-yd)*(Ez(1)*(1.d0-xd)+Ez(2)*xd)+yd*(Ez(3) *(1.d0-xd)
     $        +Ez(4)*xd)

      else


c     1 coin
         ii=imin+nx*(jmin-1)+nxy*(kmin-1)

c     WRITE(*,*) 'rr',(xs(i)-xmin)/dx,imin,jmin,kmin,ii

         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(1)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(1)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(1)=Ezr*uncomp+icomp*Ezi

c     2 coin
         
         ii=imin+1+nx*(jmin-1)+nxy*(kmin-1)
c     WRITE(*,*) 'rr2',ii
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(2)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(2)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(2)=Ezr*uncomp+icomp*Ezi

c     3 coin
         ii=imin+nx*jmin+nxy*(kmin-1)
c     WRITE(*,*) 'rr3',ii
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(3)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(3)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(3)=Ezr*uncomp+icomp*Ezi

c     4 coin
         ii=imin+1+nx*jmin+nxy*(kmin-1)
c     WRITE(*,*) 'rr4',ii
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(4)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(4)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(4)=Ezr*uncomp+icomp*Ezi

c     5 coin
         ii=imin+nx*(jmin-1)+nxy*kmin
c     WRITE(*,*) 'rr5',ii
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(5)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(5)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(5)=Ezr*uncomp+icomp*Ezi

c     6 coin
         ii=imin+1+nx*(jmin-1)+nxy*kmin
c     WRITE(*,*) 'rr6',ii
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(6)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(6)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(6)=Ezr*uncomp+icomp*Ezi

c     7 coin
         ii=imin+nx*jmin+nxy*kmin
c     WRITE(*,*) 'rr7',ii
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(7)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(7)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(7)=Ezr*uncomp+icomp*Ezi

c     8 coin
         ii=imin+1+nx*jmin+nxy*kmin
c     WRITE(*,*) 'rr8',ii
         read(11,FMT='(D22.15)',rec=ii) Exr
         read(12,FMT='(D22.15)',rec=ii) Exi
         Ex(8)=Exr*uncomp+icomp*Exi
         read(13,FMT='(D22.15)',rec=ii) Eyr
         read(14,FMT='(D22.15)',rec=ii) Eyi
         Ey(8)=Eyr*uncomp+icomp*Eyi
         read(15,FMT='(D22.15)',rec=ii) Ezr
         read(16,FMT='(D22.15)',rec=ii) Ezi
         Ez(8)=Ezr*uncomp+icomp*Ezi

c     interpolation trilineaire      
c     c_{00} = V[x_0,y_0, z_0] (1 - x_d) + V[x_1, y_0, z_0] x_d 
c     c_{10} = V[x_0,y_1, z_0] (1 - x_d) + V[x_1, y_1, z_0] x_d 
c     c_{01} = V[x_0,y_0, z_1] (1 - x_d) + V[x_1, y_0, z_1] x_d 
c     c_{11} = V[x_0,y_1, z_1] (1 - x_d) + V[x_1, y_1, z_1] x_d 
c     c_0 = c_{00}(1 - y_d) + c_{10}y_d
c     c_1 = c_{01}(1 - y_d) + c_{11}y_d
c     c = c_0(1 - z_d) + c_1z_d .    

         Exx=(1.d0-zd)*((1.d0-yd)*(Ex(1)*(1.d0-xd) +Ex(2)*xd) +yd
     $        *(Ex(3)*(1.d0-xd) + Ex(4)*xd )) +zd*((Ex(5)*(1.d0-xd) +
     $        Ex(6)*xd)*(1.d0-yd) + (Ex(7)*(1.d0-xd) + Ex(8)*xd)*yd)

         Eyy=(1.d0-zd)*((1.d0-yd)*(Ey(1)*(1.d0-xd) +Ey(2)*xd) +yd
     $        *(Ey(3)*(1.d0-xd) + Ey(4)*xd )) +zd*((Ey(5)*(1.d0-xd) +
     $        Ey(6)*xd)*(1.d0-yd) + (Ey(7)*(1.d0-xd) + Ey(8)*xd)*yd)

         Ezz=(1.d0-zd)*((1.d0-yd)*(Ez(1)*(1.d0-xd) +Ez(2)*xd) +yd
     $        *(Ez(3)*(1.d0-xd) + Ez(4)*xd )) +zd*((Ez(5)*(1.d0-xd) +
     $        Ez(6)*xd)*(1.d0-yd) + (Ez(7)*(1.d0-xd) + Ez(8)*xd)*yd)

      endif

      
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      end
