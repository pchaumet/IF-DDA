      subroutine objetpara(trope,eps,epsani,eps0,xs,ys,zs,k0,aretecube
     $     ,tabdip,nnnr,nmax,nbsphere,ndipole,nx,ny,nz,nxm,nym,nzm
     $     ,methode,epsilon ,polarisa,sidex,sidey,sidez,xg,yg,zg,phi
     $     ,theta,psi,na,nmat,file_id ,group_iddip,infostr,nstop)

#ifdef USE_HDF5
      use HDF5
#endif

      implicit none
      integer nmax,tabdip(nmax),nbsphere,ndipole,nx,ny,nz,ii,jj,i,j,k
     $     ,test,IP(3),nnnr,dddis,inv,na,nstop,nxm,nym,nzm,nmat
      double precision xs(nmax),ys(nmax),zs(nmax),k0,xg,yg,zg,x,y,z,eps0
     $     ,aretecube,sidex,sidey,sidez,side,phi,theta,psi,pi,mat(3,3)
     $     ,xr,yr,zr,x1,x2,y1,y2,z1,z2
      double complex eps,epsani(3,3),polaeps(3,3),polarisa(nmax,3,3)
     $     ,epsilon(nmax,3,3),ctmp

      character*2 methode
      character*3 trope
      character(64) infostr

      character(LEN=100) :: datasetname
#ifndef USE_HDF5
      integer,parameter:: hid_t=4
#endif

      integer(hid_t) :: file_id
      integer(hid_t) :: group_iddip
      integer :: dim(4)
      integer error
      
      write(*,*) 'object2'
c     Initialization
      nbsphere=0
      ndipole=0 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO SCHEDULE(STATIC)
      do i=1,nmax
         Tabdip(i)=0
         polarisa(i,1,1)=0.d0
         polarisa(i,1,2)=0.d0
         polarisa(i,1,3)=0.d0
         polarisa(i,2,1)=0.d0
         polarisa(i,2,2)=0.d0
         polarisa(i,2,3)=0.d0
         polarisa(i,3,1)=0.d0
         polarisa(i,3,2)=0.d0
         polarisa(i,3,3)=0.d0
         epsilon(i,1,1)=0.d0
         epsilon(i,1,2)=0.d0
         epsilon(i,1,3)=0.d0
         epsilon(i,2,1)=0.d0
         epsilon(i,2,2)=0.d0
         epsilon(i,2,3)=0.d0
         epsilon(i,3,1)=0.d0
         epsilon(i,3,2)=0.d0
         epsilon(i,3,3)=0.d0
      enddo
!$OMP ENDDO 
!$OMP END PARALLEL     
      dddis=1
      inv=1
      pi=dacos(-1.d0)

c     mesh 
      open(20,file='x.mat')
      open(21,file='y.mat')
      open(22,file='z.mat')  
c     discretization of the object under study
      open(10,file='xc.mat')
      open(11,file='yc.mat')
      open(12,file='zc.mat')  

      sidex=sidex*1.d-9
      sidey=sidey*1.d-9
      sidez=sidez*1.d-9
      xg=xg*1.d-9
      yg=yg*1.d-9
      zg=zg*1.d-9
      theta=theta*pi/180.d0
      phi=phi*pi/180.d0
      psi=psi*pi/180.d0

      write(*,*) 'object4'
      if (sidex.eq.0.d0) then
         nstop=1
         infostr='object cuboid: sidex=0!'
         return
      elseif (sidey.eq.0.d0) then
         nstop=1
         infostr='object cuboid: sidey=0!'
         return
      elseif (sidez.eq.0.d0) then
         nstop=1
         infostr='object cuboid: sidez=0!'
         return
      endif
      write(*,*) 'object'


      mat(1,1)=dcos(psi)*dcos(phi)-dsin(psi)*dcos(theta)*dsin(phi)
      mat(1,2)=-dcos(psi)*dsin(phi)-dsin(psi)*dcos(theta)*dcos(phi)
      mat(1,3)=dsin(psi)*dsin(theta)
      
      mat(2,1)=dsin(psi)*dcos(phi)+dcos(psi)*dcos(theta)*dsin(phi)
      mat(2,2)=-dsin(psi)*dsin(phi)+dcos(psi)*dcos(theta)*dcos(phi)
      mat(2,3)=-dcos(psi)*dsin(theta)

      mat(3,1)=dsin(theta)*dsin(phi)
      mat(3,2)=dsin(theta)*dcos(phi)
      mat(3,3)=dcos(theta)

      x1=1.d30
      x2=-1.d30
      y1=1.d30
      y2=-1.d30
      z1=1.d30
      z2=-1.d30

      x=-sidex/2.d0
      y=-sidey/2.d0
      z=-sidez/2.d0
      xr=mat(1,1)*x+mat(1,2)*y+mat(1,3)*z
      yr=mat(2,1)*x+mat(2,2)*y+mat(2,3)*z
      zr=mat(3,1)*x+mat(3,2)*y+mat(3,3)*z
      x1=min(x1,xr)
      x2=max(x2,xr)
      y1=min(y1,yr)
      y2=max(y2,yr)
      z1=min(z1,zr)
      z2=max(z2,zr)

      x=sidex/2.d0
      y=-sidey/2.d0
      z=-sidez/2.d0
      xr=mat(1,1)*x+mat(1,2)*y+mat(1,3)*z
      yr=mat(2,1)*x+mat(2,2)*y+mat(2,3)*z
      zr=mat(3,1)*x+mat(3,2)*y+mat(3,3)*z
      x1=min(x1,xr)
      x2=max(x2,xr)
      y1=min(y1,yr)
      y2=max(y2,yr)
      z1=min(z1,zr)
      z2=max(z2,zr)

      x=-sidex/2.d0
      y=sidey/2.d0
      z=-sidez/2.d0
      xr=mat(1,1)*x+mat(1,2)*y+mat(1,3)*z
      yr=mat(2,1)*x+mat(2,2)*y+mat(2,3)*z
      zr=mat(3,1)*x+mat(3,2)*y+mat(3,3)*z
      x1=min(x1,xr)
      x2=max(x2,xr)
      y1=min(y1,yr)
      y2=max(y2,yr)
      z1=min(z1,zr)
      z2=max(z2,zr)

      x=-sidex/2.d0
      y=-sidey/2.d0
      z=sidez/2.d0
      xr=mat(1,1)*x+mat(1,2)*y+mat(1,3)*z
      yr=mat(2,1)*x+mat(2,2)*y+mat(2,3)*z
      zr=mat(3,1)*x+mat(3,2)*y+mat(3,3)*z
      x1=min(x1,xr)
      x2=max(x2,xr)
      y1=min(y1,yr)
      y2=max(y2,yr)
      z1=min(z1,zr)
      z2=max(z2,zr)

      x=sidex/2.d0
      y=sidey/2.d0
      z=-sidez/2.d0
      xr=mat(1,1)*x+mat(1,2)*y+mat(1,3)*z
      yr=mat(2,1)*x+mat(2,2)*y+mat(2,3)*z
      zr=mat(3,1)*x+mat(3,2)*y+mat(3,3)*z
      x1=min(x1,xr)
      x2=max(x2,xr)
      y1=min(y1,yr)
      y2=max(y2,yr)
      z1=min(z1,zr)
      z2=max(z2,zr)


      x=sidex/2.d0
      y=-sidey/2.d0
      z=sidez/2.d0
      xr=mat(1,1)*x+mat(1,2)*y+mat(1,3)*z
      yr=mat(2,1)*x+mat(2,2)*y+mat(2,3)*z
      zr=mat(3,1)*x+mat(3,2)*y+mat(3,3)*z
      x1=min(x1,xr)
      x2=max(x2,xr)
      y1=min(y1,yr)
      y2=max(y2,yr)
      z1=min(z1,zr)
      z2=max(z2,zr)

      x=-sidex/2.d0
      y=sidey/2.d0
      z=sidez/2.d0
      xr=mat(1,1)*x+mat(1,2)*y+mat(1,3)*z
      yr=mat(2,1)*x+mat(2,2)*y+mat(2,3)*z
      zr=mat(3,1)*x+mat(3,2)*y+mat(3,3)*z
      x1=min(x1,xr)
      x2=max(x2,xr)
      y1=min(y1,yr)
      y2=max(y2,yr)
      z1=min(z1,zr)
      z2=max(z2,zr)

      x=sidex/2.d0
      y=sidey/2.d0
      z=sidez/2.d0
      xr=mat(1,1)*x+mat(1,2)*y+mat(1,3)*z
      yr=mat(2,1)*x+mat(2,2)*y+mat(2,3)*z
      zr=mat(3,1)*x+mat(3,2)*y+mat(3,3)*z

      x1=min(x1,xr)
      x2=max(x2,xr)
      y1=min(y1,yr)
      y2=max(y2,yr)
      z1=min(z1,zr)
      z2=max(z2,zr)

      side=max(x2-x1,y2-y1,z2-z1)
      aretecube=side/dble(nnnr)


      call inversemat33r(mat)

c     boite au plus pres de l'objet
      nx=int((x2-x1)/aretecube)
      ny=int((y2-y1)/aretecube)
      nz=int((z2-z1)/aretecube)
      write(*,*) 'nnn',nx,ny,nz

      if (nx.gt.nxm.or.ny.gt.nym.or.nz.gt.nzm) then
         nstop=1
         infostr='Dimension Problem of the Box : Box too small!'
         write(99,*) 'dimension Problem',nx,nxm,ny,nym,nz,nzm
         write(*,*) 'dimension Problem',nx,nxm,ny,nym,nz,nzm
         return
      endif
      
      if (na.eq.-1) then

         do i=1,nz
            do j=1,ny
               do k=1,nx

                  x=-(x2-x1)/2.d0+aretecube*(dble(k)-0.5d0)
                  y=-(y2-y1)/2.d0+aretecube*(dble(j)-0.5d0)
                  z=-(z2-z1)/2.d0+aretecube*(dble(i)-0.5d0)                 

                  if (j.eq.1.and.k.eq.1.and.nmat.ne.1) write(22,*) z+zg
                  if (i.eq.1.and.k.eq.1.and.nmat.ne.1) write(21,*) y+yg
                  if (j.eq.1.and.i.eq.1.and.nmat.ne.1) write(20,*) x+xg

                  ndipole=ndipole+1                   

                  xr=mat(1,1)*x+mat(1,2)*y+mat(1,3)*z
                  yr=mat(2,1)*x+mat(2,2)*y+mat(2,3)*z
                  zr=mat(3,1)*x+mat(3,2)*y+mat(3,3)*z


                  if (dabs(xr).le.sidex/2.d0.and.dabs(yr).le.sidey
     $                 /2.d0.and.dabs(zr).le.sidez/2.d0) then  

                     nbsphere=nbsphere+1

                     Tabdip(ndipole)=nbsphere
                     xs(nbsphere)=x+xg
                     ys(nbsphere)=y+yg
                     zs(nbsphere)=z+zg                   
                     
                     if (trope.eq.'iso') then
                        call poladiff(aretecube,eps,eps0,k0,dddis
     $                       ,methode,ctmp)  
                        polarisa(nbsphere,1,1)=ctmp
                        polarisa(nbsphere,2,2)=ctmp
                        polarisa(nbsphere,3,3)=ctmp
                        epsilon(nbsphere,1,1)=eps
                        epsilon(nbsphere,2,2)=eps
                        epsilon(nbsphere,3,3)=eps
                     else
                        call polaepstens(aretecube,epsani,eps0,k0,dddis
     $                       ,methode,inv,polaeps)
                        do ii=1,3
                           do jj=1,3
                              epsilon(nbsphere,ii,jj)=epsani(ii,jj)
                              polarisa(nbsphere,ii,jj)=polaeps(ii,jj)
                           enddo
                        enddo
                     endif
                  endif
               enddo
            enddo
         enddo 
      elseif (na.eq.0) then
         

c     fin essai
         do i=1,nz
            do j=1,ny
               do k=1,nx
c     x=-side/2.d0+aretecube*(dble(k)-0.5d0)
c     y=-side/2.d0+aretecube*(dble(j)-0.5d0)
c     z=-side/2.d0+aretecube*(dble(i)-0.5d0)
                  x=-(x2-x1)/2.d0+aretecube*(dble(k)-0.5d0)
                  y=-(y2-y1)/2.d0+aretecube*(dble(j)-0.5d0)
                  z=-(z2-z1)/2.d0+aretecube*(dble(i)-0.5d0)                 

                  

                  if (j.eq.1.and.k.eq.1.and.nmat.ne.1) write(22,*) z+zg
                  if (i.eq.1.and.k.eq.1.and.nmat.ne.1) write(21,*) y+yg
                  if (j.eq.1.and.i.eq.1.and.nmat.ne.1) write(20,*) x+xg

                  ndipole=ndipole+1  
                  nbsphere=nbsphere+1
                  Tabdip(ndipole)=nbsphere
                  
                  xr=mat(1,1)*x+mat(1,2)*y+mat(1,3)*z
                  yr=mat(2,1)*x+mat(2,2)*y+mat(2,3)*z
                  zr=mat(3,1)*x+mat(3,2)*y+mat(3,3)*z

                  if (dabs(xr).le.sidex/2.d0.and.dabs(yr).le.sidey
     $                 /2.d0.and.dabs(zr).le.sidez/2.d0) then
     $                 
                     xs(nbsphere)=x+xg
                     ys(nbsphere)=y+yg
                     zs(nbsphere)=z+zg
                     if (trope.eq.'iso') then
                        call poladiff(aretecube,eps,eps0,k0,dddis
     $                       ,methode,ctmp)  
                        polarisa(nbsphere,1,1)=ctmp
                        polarisa(nbsphere,2,2)=ctmp
                        polarisa(nbsphere,3,3)=ctmp
                        epsilon(nbsphere,1,1)=eps
                        epsilon(nbsphere,2,2)=eps
                        epsilon(nbsphere,3,3)=eps

                     else
                        call polaepstens(aretecube,epsani,eps0,k0,dddis
     $                       ,methode,inv,polaeps)
                        do ii=1,3
                           do jj=1,3
                              epsilon(nbsphere,ii,jj)=epsani(ii,jj)
                              polarisa(nbsphere,ii,jj)=polaeps(ii,jj)
                           enddo
                        enddo
                     endif
                  else
                     xs(nbsphere)=x+xg
                     ys(nbsphere)=y+yg
                     zs(nbsphere)=z+zg     
                     epsilon(nbsphere,1,1)=(1.d0,0.d0)
                     epsilon(nbsphere,2,2)=(1.d0,0.d0)
                     epsilon(nbsphere,3,3)=(1.d0,0.d0)
                  endif
               enddo
            enddo
         enddo
      else 
         infostr='na should be equal to -1 or 0'
         nstop=1
         return
      endif

      if (ndipole.gt.nmax) then
         infostr='nmax parameter too small: increase nxm nym nzm'
         nstop=1
         return
      endif

      if (nmat.eq.0) then
         do i=1,nbsphere
            write(10,*) xs(i)
            write(11,*) ys(i)
            write(12,*) zs(i)
         enddo
      elseif (nmat.eq.2) then
         dim(1)=nbsphere
         dim(2)=nmax
         datasetname="Dipole position x"
         call hdf5write1d(group_iddip,datasetname,xs,dim)
         datasetname="Dipole position y"
         call hdf5write1d(group_iddip,datasetname,ys,dim)
         datasetname="Dipole position z"
         call hdf5write1d(group_iddip,datasetname,zs,dim)
      endif
      
      close(10)
      close(11)
      close(12)
      close(20)
      close(21)
      close(22)

      write(*,*) 'Cuboid::nbsphere=',nbsphere
      end
